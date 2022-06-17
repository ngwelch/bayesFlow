#!/usr/bin/env Rscript

startTime = Sys.time()
tm = as.character(Sys.time())
cat("\nStarting MCMC Sampler", tm, "\n" )


library(data.table)
library(magrittr)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(MASS)
library(nimble)

args = commandArgs(TRUE) # args=c(50, 0.5, 61)

dir = as.character( args[1] )
origin.number = as.integer( args[2] )

#ss.mcmc = as.integer( args[2] )
#thin.mcmc = as.integer( args[3] )
#trace.size = as.integer( args[4] )

# use ss.mcmc and thin.mcmc parameters as the outflow model
load(file=paste0(dir, "/src/mcmc_params.RData") )


## --------------------------------------------------------------------------------------
tm = as.character(Sys.time())
cat("\nLoading data", tm , "\n" )

unlargest200 = fread(file=paste0(dir, "/data/200isoRegionCodes.csv"))[order(iso)]
isos = unlargest200[order(iso)][,iso]

#saveRDS(trace.size, file=paste0(dir, "/heldout/traceSize.rds") )

thisOrigin = isos[origin.number]
thisOrigCode = unlargest200[iso==thisOrigin, country_code]

unlink(paste0(dir, "/heldout/kappa/", thisOrigCode,".csv"))
unlink(paste0(dir, "/heldout/psi/", thisOrigCode,".csv"))

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



############################# ORIGIN DATA ########################################

flowTrainOrigin = flowTrain[orig==thisOrigin][,.(dest, flow, year0)]
tmp = dcast(flowTrainOrigin, formula=dest~year0, value.var="flow")

m = as.matrix( tmp[,-1]) %>% set_rownames( tmp$dest )
N = colSums( m )

nt = ncol( m )
nd = nrow( m )

destIsos = flowTrainOrigin[,unique(dest)] %>% sort()

############################### BHM MODEL ###########################################

tm = as.character(Sys.time())
cat("\nInitializing model", tm, "\n" )


flowCode <- nimbleCode({
  # likelihood #
  for( t in 1:Nt ){
  	m[1:J, t] ~ dmulti( size=N[t], prob=p[1:J,t] )
  } # end t
  # destination component #
  for( t in 1:Nt ){
	for( j in 1:J ){
		eta[j, t] ~ dnorm( mean=kappa[j], sd=psi[j] )
		p[j,t] <- exp( eta[j,t] ) / sum( exp( eta[1:J,t] ) )
    } # end t
  } # end j
  # hyperpriors
  for(j in 1:J){
	kappa[j] ~ dnorm(mean=0, sd=10)
    psi[j] ~ dbeta(p0, q0)
  } # end j
})

p0 = 7.409950 #5.755366 #4.404856 #3.44402
q0 = 3.893723 #2.678973 #3.121439 #6.61646

pi0 = apply( m+0.25, 2, function(x) x/sum(x) )
eta0 = apply( m+ 0.25, 2, function(x) log( x / exp( 1/length(x) * sum( log(x) ) ) ) )
kappa0.tmp = rowMeans( eta0 )
kappa0 = kappa0.tmp - mean(kappa0.tmp)
psi0 = rbeta(n=nd, p0, q0)

flow <- nimbleModel(code = flowCode,
                    name = "flow",
                    constants = list( J=nd, Nt=nt, p0=p0, q0=q0),
                    data = list( N=N,m=m ),
                    inits = list( kappa=kappa0, psi=psi0, eta=eta0, p=pi0 ),
                    calculate=TRUE )

flowConf <- configureMCMC(model = flow,
                          monitors = c("kappa", "psi"),
                          thin = thin.mcmc,
                          useConjugacy = TRUE,
                          print = FALSE)

Cflow <- compileNimble( flow )
flowMCMC <- buildMCMC( flowConf )
CflowMCMC <- compileNimble(flowMCMC, project=flow, resetFunctions=FALSE )

################################# BHM SAMPLE #########################################

tm = as.character(Sys.time())
cat("\nStarting sampler at", tm, "\n" )

CflowMCMC$run( 100 )
CflowMCMC$run( ss.mcmc )
postSample = CflowMCMC$mvSamples %>% as.matrix() %>% as.data.table()
thin.ss = nrow( postSample )

isoj = data.table(j=1:nd, dest=destIsos)

# kappa
kappa.mcmc = 
  postSample[,.SD, .SDcols=which( str_detect(colnames(postSample), "kappa\\["))] %>%
  t() %>%
  as.data.table(keep.rownames=TRUE) %>%
  set_colnames(c("j", 1:thin.ss)) %>%
  .[,j:=tstrsplit(j, "\\[", keep=2)] %>%
  .[,j:=tstrsplit(j, "\\]", keep=1)] %>%
  .[,j:=as.numeric(j)] %>%
  merge(isoj, by="j") %>%
  setcolorder(c("dest", "j"))
fwrite(kappa.mcmc, file=paste0(dir, "/heldout/kappa/", thisOrigCode ,".csv"))

# psi 
psi.mcmc = 
  postSample[,.SD, .SDcols=which( str_detect(colnames(postSample), "psi\\["))] %>%
  t() %>%
  as.data.table(keep.rownames=TRUE) %>%
  set_colnames(c("j", 1:thin.ss)) %>%
  .[,j:=tstrsplit(j, "\\[", keep=2)] %>%
  .[,j:=tstrsplit(j, "\\]", keep=1)] %>%
  .[,j:=as.numeric(j)] %>%
  merge(isoj, by="j") %>%
  setcolorder(c("dest", "j"))
fwrite(psi.mcmc, file=paste0(dir, "/heldout/psi/", thisOrigCode,".csv"))


tm = as.character(Sys.time())
cat("\n", thisOrigin, "complete", tm, "\n" )

