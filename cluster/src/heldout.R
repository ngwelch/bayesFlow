#!/usr/bin/env Rscript


startTime = Sys.time()
tm = as.character(Sys.time())
cat("\nStarting sampler evaluation", tm, "\n" )


library(data.table)
library(magrittr)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(nimble)
theme_set(theme_cowplot())
library(MASS)

args = commandArgs(TRUE) # args=c("~/research/raftery/welchRaftery2022", "0.5")

dir = as.character( args[1] )
burnin.proportion = as.numeric( args[2] )
#trace.annual = round( readRDS(file=paste0(dir, "/heldout/traceSize.rds") ) / 5 )

## --------------------------------------------------------------------------------------

unlargest200 = fread(paste0(dir, "/data/200isoRegionCodes.csv"))[order(iso)]
isos = unlargest200[order(iso)][,iso] #c("DEU", "FRA", "GBR", "IRL", "NLD") 

tm = as.character(Sys.time())
cat("\nLoading data at", tm , "\n" )

flowTrain = 
  fread(
    file=paste0(dir, "/data/azoseRaftery2019flows.csv")
    )[,.(orig_code=origin, 
         dest_code=destination, 
         orig=origIso, 
         dest=destIso, 
         year0=year,
         flow=migrantCount)][
           orig %in% isos & dest %in% isos]

# marginal flows
marginflowTrain = 
  merge(flowTrain[,.(inflow=sum(flow)), .(year0, iso=dest)], 
        flowTrain[,.(outflow=sum(flow)), .(year0, iso=orig)], 
        by=c("year0", "iso"))[,netflow:=inflow-outflow]

periodTrain = flowTrain[, sort( unique( year0 ) ) ]

nt = length( periodTrain )
nc = length( isos )
nd = nc - 1

## --------------------------------------------------------------------------------------
# Population
data("pop", package="wpp2015")
data("popproj", package="wpp2015")

wppHistoricPop = as.data.table(pop)
wppProjectedPop = as.data.table(popproj)

popdt = 
  merge(wppHistoricPop, wppProjectedPop, by=c("name", "country_code")) %>%
  merge(unlargest200, by="country_code")  %>%
  setnames(old="name", new="country_name") %>%
  setcolorder(c("area_code", "reg_code", "country_code", 
                "area_name", "reg_name", "country_name", "iso"))

popNames = c("iso", as.character(seq(1950, 2015, 5)))
popIsodt = 
  data.table::melt(popdt[,eval(popNames), with=FALSE], 
       id.vars="iso", variable.name="year0", value.name="midPeriodPop000", 
       variable.factor = FALSE)[,year0:=as.numeric(year0)][iso %in% isos][
         ,`:=`(pop=midPeriodPop000*1e3, midPeriodPop000=NULL)
       ][]



## --------------------------------------------------------------------------------------
popflowTrain = merge(marginflowTrain, popIsodt)[
  ,`:=`(netrate=netflow/(5*pop/1e3), 
        outrate=outflow/(5*pop/1e3))]

tmp = 
  popflowTrain[,.(year0, iso, outrate)] %>%
  data.table::dcast(iso ~ year0, value.var="outrate", fill=0)

outrateTrain = 
  tmp[,-"iso"] %>%
  as.matrix() %>%
  set_rownames(tmp$iso)


tmp = marginflowTrain[,.(iso, year0, outflow)] %>% dcast(formula=iso~year0)

outflowTrain = 
  tmp[,-"iso"] %>%
  as.matrix() %>%
  set_rownames(tmp$iso)


## --------------------------------------------------------------------------------------
data("pop", package="wpp2019")

pop2015 = 
  merge(unique( flowTrain[,.(country_code=orig_code, iso=orig)] ), 
        data.table(country_code=pop[,"country_code"], 
                   pop=1e3*pop[,"2015"], 
                   year0=2015) ) %>%
  .[,country_code:=NULL] %>%
  .[]


flowTest = 
  fread(file=paste0(dir, "/data/abelCohen2019flowsv6_flowdt.csv") )[year0==2015][orig %in% isos & dest %in% isos]
marginflowTest = 
  merge(flowTest[,.(inflow=sum(flow)), .(year0, iso=dest)], 
        flowTest[,.(outflow=sum(flow)), .(year0, iso=orig)], 
        by=c("year0", "iso"))[,netflow:=inflow-outflow]


popflowTest = merge(marginflowTest, pop2015)[
  ,`:=`(netrate=netflow/(5*pop/1e3), 
        outrate=outflow/(5*pop/1e3))]

outrateTest = popflowTest[,.(year0, iso, outrate)]
piTest = flowTest[,.(dest, persist=flow/sum(flow)),orig][order(orig, dest)]


############################ OUT OF SAMPLE EVALUATION #################################

## pi predictions ##
tm = as.character(Sys.time())
cat("\nGenerating inflow predictions", tm, "\n")

# MCMC outflow parameter sample #
mu.all = fread( file=file.path(dir, "heldout", "mu.csv"), header=TRUE )
sigma.all = fread( file=file.path(dir, "heldout", "sigma.csv"), header=TRUE )

postSampleSize = ncol( mu.all ) - 2 
postSampleIndex = floor( burnin.proportion * ( postSampleSize ) + 1 ) : postSampleSize
outflowNames = c( "orig", as.character( postSampleIndex ) )
inflowNames = c( "dest", as.character( postSampleIndex ) )

S = length(postSampleIndex)  
mcmcNames = c("orig", as.character( 1:S ) )

mu.post = mu.all[,..outflowNames] %>% set_colnames( mcmcNames )
sigma.post = sigma.all[,..outflowNames] %>% set_colnames( mcmcNames )
phi.post = fread( file=file.path(dir, "heldout", "phi.csv"), header=TRUE )[postSampleIndex,phi]


piPredict = 
  array(0, 
        dim=c(nc, nc, S), 
        dimnames=list(orig=isos, dest=isos, sample=1:S)
  )


for(o in isos){
  cat("\n", o )
  occ = unlargest200[iso==o, country_code]
  kappa.post = fread( file=paste0(dir, "/heldout/kappa/", occ,".csv"), header=TRUE )[,..inflowNames]
  psi.post = fread( file=paste0(dir, "/heldout/psi/", occ,".csv"), header=TRUE )[,..inflowNames]
  destIsos = kappa.post[,dest] 
  nd = length( destIsos )
  eta = matrix(0, nrow=S, ncol=nd, dimnames=list(1:S, destIsos))
  for(d in destIsos){
    kappa.sample = as.numeric( kappa.post[dest==d][,-"dest"] )
    psi.sample = as.numeric( psi.post[dest==d][,-"dest"] )
    eta[,d] = rnorm(n=S, mean=kappa.sample, sd=psi.sample)
  }
  #eta = sweep( eta, 1, rowMeans(etaRaw) )
  
  piPredict[orig=o, destIsos,] = apply(eta, 1, function(s) exp(s)/sum( exp(s) ) )
}
# use Sudan to approximate South Sudan flow preferences
piPredict[orig="SSD",,] = piPredict[orig="SDN",,]
piPredict[orig="SSD", dest="SDN",] = piPredict[orig="SDN", dest="SSD",] 
piPredict[orig="SSD", dest="SSD",] = 0


## flow predictions ##
tm = as.character(Sys.time())
cat("\nGenerating flow predictions", tm, "\n")

flowPredict = destSharePredict = 
  array(0, 
        dim=c(nc, nc, S),
        dimnames=list(origin=isos, destination=isos, sample=1:S)
  )

outflowPredict = matrix(nrow=nc, ncol=S, dimnames=list( isos, 1:S) )

for(o in isos){
  cat("\n", o )
  logDeltaLast = log( outrateTrain[o,"2010"] )
  pop.orig = pop2015[iso==o, pop]
  
  mu.sample.o = as.numeric( mu.post[orig==o][,-"orig"] )
  sigma.sample.o = as.numeric( sigma.post[orig==o][,-"orig"] )
  rho.sample.o = (1-phi.post)*mu.sample.o + phi.post*logDeltaLast
  
  logDeltaNext = rnorm(n=S, mean=rho.sample.o, sd=sigma.sample.o)
  deltaNext = exp( logDeltaNext )
  
  outflowPredict[o,] = ceiling( deltaNext * popdt[iso==o, `2015`] * 5 )
  
  destIsos = isos[isos!=o]
  
  for(s in 1:S){
    outflowPredict.s = outflowPredict[o,s]
    piPredict.s = piPredict[orig=o, dest=destIsos, s]
    flowPredict[o,destIsos,s] = nimble::rmulti(n=1, size=outflowPredict.s, prob=piPredict.s)
  }
}

## pit estimate ##
tm = as.character(Sys.time())
cat("\nGenerating PIT results", tm, "\n")


pit.flow = array(NA, dim=c(nc, nc), dimnames=list(orig=isos, dest=isos))

for(o in isos){
  for(d in isos){
    if(o==d) next
    mig.true = flowTest[orig==o & dest==d, flow]
    mig.sample = flowPredict[o,d,]
	
	p0 = ecdf( mig.sample )( mig.true - 1 )
	p1 = ecdf( mig.sample )( mig.true )
	v = runif(1)
	u = p0 + v*(p1-p0)
	pit.flow[o,d] = u
    
  }
}
pit.flow = reshape2::melt(pit.flow, value.name="pit")
pit.flow = pit.flow[complete.cases(pit.flow),]

flowCoverage = 
  as.data.table(pit.flow)[,.("80%"=round(100*mean(pit>=0.1 & pit<=0.9)),
                             "90%"=round(100*mean(pit>=0.05 & pit<=0.95)),
                             "95%"=round(100*mean(pit>=0.025 & pit<=0.975))
  )]


# outflows 
pit.outflow = rep(NA, nc) %>% set_names(isos)
mig.outflow.sample = apply(flowPredict, MARGIN=c("sample", "origin"), sum) 
for(o in isos){
  mig.true = marginflowTest[iso==o, outflow]
  mig.sample = mig.outflow.sample[,o]
  
  p0 = ecdf( mig.sample )( mig.true - 1 )
  p1 = ecdf( mig.sample )( mig.true )
  v = runif(1)
  u = p0 + v*(p1-p0)
  
  pit.outflow[o] = u
}

outflowCoverage = 
  data.table("80%"=round(100*mean(pit.outflow>=0.1 & pit.outflow<=0.9)),
             "90%"=round(100*mean(pit.outflow>=0.05 & pit.outflow<=0.95)),
             "95%"=round(100*mean(pit.outflow>=0.025 & pit.outflow<=0.975))
              )

# inflows
pit.inflow = rep(NA, nc) %>% set_names( isos )
mig.inflow.sample = apply(flowPredict, MARGIN=c("sample", "destination"), sum) 
for(d in isos){
  mig.true = marginflowTest[iso==d, inflow]
  mig.sample = mig.inflow.sample[,d]
  
  p0 = ecdf( mig.sample )( mig.true - 1 )
  p1 = ecdf( mig.sample )( mig.true )
  v = runif(1)
  u = p0 + v*(p1-p0)
  
  pit.inflow[d] = u
}

inflowCoverage = 
  data.table("80%"=round(100*mean(pit.inflow>=0.1 & pit.inflow<=0.9)),
             "90%"=round(100*mean(pit.inflow>=0.05 & pit.inflow<=0.95)),
             "95%"=round(100*mean(pit.inflow>=0.025 & pit.inflow<=0.975))
  )


# netflows 
pit.netflow = rep(NA, nc) %>% set_names(isos)
mig.netflow.sample = mig.inflow.sample - mig.outflow.sample

for(thisIso in isos){
  mig.true = marginflowTest[iso==thisIso, netflow]
  mig.sample = mig.netflow.sample[,thisIso] 
  
  p0 = ecdf( mig.sample )( mig.true - 1 )
  p1 = ecdf( mig.sample )( mig.true )
  v = runif(1)
  u = p0 + v*(p1-p0)
  
  pit.netflow[thisIso] = u
}

netflowCoverage = 
  data.table("80%"=round(100*mean(pit.netflow>=0.1 & pit.netflow<=0.9)),
             "90%"=round(100*mean(pit.netflow>=0.05 & pit.netflow<=0.95)),
             "95%"=round(100*mean(pit.netflow>=0.025 & pit.netflow<=0.975))
  )

coverageSummary = 
  rbind("flow"=flowCoverage, 
        "out"=outflowCoverage, 
        "in"=inflowCoverage,
        "net"=netflowCoverage) %>%
  cbind( data.table( "Margin"=c("flow", "out", "in", "net") ) ) %>%
  setcolorder("Margin")

pit.margin = cbind( "out"=pit.outflow, "in"=pit.inflow, "net"=pit.netflow ) 

## predictions ##

outflowPost = as.data.table( mig.outflow.sample )
outflowMedian = apply( outflowPost, 2, median )
outflowObs = marginflowTest[,.(iso, obs=outflow)]

inflowPost = as.data.table( mig.inflow.sample )
inflowMedian = apply( inflowPost, 2, median )
inflowObs = marginflowTest[,.(iso, obs=inflow)]

netflowPost = as.data.table( mig.netflow.sample )
netflowMedian = apply( netflowPost, 2, median )
netflowObs = marginflowTest[,.(iso, obs=netflow)]

## error 
tm = as.character(Sys.time())
cat("\nCompiling error results", tm, "\n")


# flow error
flowPersist = flowTrain[year0==2010,.(orig, dest, persist=flow)]
flowHistMean = flowTrain[,.(histMean=round(mean(flow))), .(orig, dest)]
flowBayes = 
  apply(flowPredict, MARGIN=c("origin", "destination"), function(x) round(median(x)) ) %>% 
  reshape2::melt(value.name = "bayes") %>% 
  as.data.table() %>%
  .[origin!=destination] %>%
  setnames(old=c("origin", "destination"), new=c("orig", "dest"))

flowMerge = 
  Reduce(merge, 
         list(flowPersist, flowHistMean, flowBayes, flowTest[,.(orig, dest, obs=flow)]))

flowMAE = 
  flowMerge[,.(persist=mean(abs(persist-obs)), 
               histMean=mean(abs(histMean-obs)),
               bayes=mean(abs(bayes-obs))
  )] %>% round() 

# outflow error
outflowPersist = flowPersist[,.(persistOut=sum(persist)),.(iso=orig)]
outflowHistMean = flowHistMean[,.(histMeanOut=sum(histMean)),.(iso=orig)]
outflowBayes = 
  apply( outflowPost, 2, function(x) round(median(x)) ) %>% 
  as.data.table(TRUE) %>%
  set_colnames(c("iso", "bayesMedianOut"))

outflowMerge = 
  Reduce(merge, 
         list(outflowPersist, outflowHistMean, outflowBayes, outflowObs))

outflowMAE = 
  outflowMerge[,.(persist=mean(abs(persistOut-obs)), 
                  histMean=mean(abs(histMeanOut-obs)),
                  bayes=mean(abs(bayesMedianOut-obs))
  )] %>% round() 

# inflow error
inflowPersist = flowPersist[,.(persistIn=sum(persist)),.(iso=dest)]
inflowHistMean = flowHistMean[,.(histMeanIn=sum(histMean)),.(iso=dest)]
inflowBayes = 
  apply( inflowPost, 2, function(x) round(median(x)) ) %>% 
  as.data.table(TRUE) %>%
  set_colnames(c("iso", "bayesMedianIn"))

inflowMerge = 
  Reduce(merge, 
         list(inflowPersist, inflowHistMean, inflowBayes, inflowObs))

inflowMAE = 
  inflowMerge[,.(persist=mean(abs(persistIn-obs)), 
                 histMean=mean(abs(histMeanIn-obs)),
                 bayes=mean(abs(bayesMedianIn-obs))
  )] %>% round() 


# net flow error
netflowPersist = 
  merge( inflowPersist, outflowPersist )[,.(iso, persistNet=persistIn-persistOut)]
netflowHistMean = 
  merge( inflowHistMean, outflowHistMean )[,.(iso, histMeanNet=histMeanIn-histMeanOut)]
netflowBayes = 
  apply( netflowPost, 2, function(x) round(median(x)) ) %>% 
  as.data.table(TRUE) %>%
  set_colnames(c("iso", "bayesMedianNet"))

netflowMerge = 
  Reduce(merge, 
         list(netflowPersist, netflowHistMean, netflowBayes, netflowObs))

netflowMAE = 
  netflowMerge[,.(persist=mean(abs(persistNet-obs)), 
                  histMean=mean(abs(histMeanNet-obs)),
                  bayes=mean(abs(bayesMedianNet-obs))
  )] %>% round() 


MAE = 
  cbind(margin=c("flow", "out", "in", "net"), 
        rbind(flows=flowMAE, 
              outflows=outflowMAE, 
              inflows=inflowMAE, 
              netflows=netflowMAE) ) %>%
  set_colnames(c("Margin", "Persistence", "Historic Mean", "BHM Median"))

## --------------------------------------------------------------------------------------

tm = as.character(Sys.time())
cat("\nWriting results", tm, "\n" )


# flow predict
saveRDS( flowPredict, file=paste0(dir, "/heldout/results/flowPredict.rds") )


# pit results
save( list=c("pit.flow", "pit.outflow", "pit.inflow", "pit.netflow"), 
      file=paste0(dir, "/heldout/results/pit.RData") )

fwrite(coverageSummary, file=paste0(dir, "/heldout/results/coverage_summary.csv") )

# MAE results
fwrite(MAE, file=paste0(dir, "/heldout/results/mae.csv") )


## --------------------------------------------------------------------------------------
tm = as.character(Sys.time())
cat("\nGenerating plots", tm, "\n" )

# pit plot
pdf( paste0(dir, "/heldout/plot/pit.pdf") )
# all pit plots 
breaks=10
hist(pit.flow$pit, breaks=seq(0, 1, 1/breaks), main="Flows", xlab="PIT")
abline(h=nc*(nc-1)/breaks)

hist(pit.outflow, breaks=seq(0, 1, 1/breaks), main="Outflows", xlab="PIT")
abline(h=nc/breaks)

hist(pit.inflow, breaks=seq(0, 1, 1/breaks), main="Inflows", xlab="PIT")
abline(h=nc/breaks)

hist(pit.netflow, breaks=seq(0, 1, 1/breaks), main="Net Flows", xlab="PIT")
abline(h=nc/breaks)
dev.off()

# marginal flows

pdf( paste0(dir, "/heldout/plot/marginal_flows.pdf") )
for(thisIso in isos){
  outflowPlot = 
    ggplot() + 
    geom_histogram(data=outflowPost, aes(x=get(thisIso), y=..density..), bins=30) + 
    geom_point(aes(x=outflowMedian[thisIso], y=0, col="Poster Median"), size=3) + 
    geom_point(data=outflowObs[iso==thisIso], 
               aes(x=obs, y=0, col="2015-2020 Observed"), size=3) +
    xlab("Outflow") + 
    ggtitle(paste(thisIso, "Marginal Flows")) +
    theme(legend.title = element_blank(), 
          legend.position="none", 
          legend.justification = "center")
  
  inflowPlot = 
    ggplot() + 
    geom_histogram(data=inflowPost, aes(x=get(thisIso), y=..density..), bins=30) + 
    geom_point(aes(x=inflowMedian[thisIso], y=0, col="Poster Median"), size=3) + 
    geom_point(data=inflowObs[iso==thisIso], 
               aes(x=obs, y=0, col="2015-2020 Observed"), size=3) +
    xlab("Inflow") + 
    theme(legend.title = element_blank(), 
          legend.position="none", 
          legend.justification = "center")
  
  netflowPlot = 
    ggplot() + 
    geom_histogram(data=netflowPost, aes(x=get(thisIso), y=..density..), bins=30) + 
    geom_point(aes(x=netflowMedian[thisIso], y=0, col="Poster Median"), size=3) + 
    geom_point(data=netflowObs[iso==thisIso], 
               aes(x=obs, y=0, col="2015-2020 Observed"), size=3) + 
    xlab("Net Flow") + 
    theme(legend.title = element_blank(), 
          legend.position="bottom", 
          legend.justification = "center")
  
  yrange = range( 0, 
                  layer_scales(outflowPlot)$y$get_limits(),
                  layer_scales(inflowPlot)$y$get_limits(),
                  layer_scales(netflowPlot)$y$get_limits() )
  
  out = plot_grid(outflowPlot + coord_cartesian(ylim=yrange), 
                  inflowPlot + coord_cartesian(ylim=yrange), 
                  netflowPlot + coord_cartesian(ylim=yrange), 
                  nrow=3, ncol=1,
                  rel_heights=c(1,0.9,1))
  print( out )
}
dev.off()

tm = as.character(Sys.time())
cat("\nCompleted heldout analysis at", tm, "\n" )



