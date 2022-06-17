#!/usr/bin/env Rscript

startTime = Sys.time()
tm = as.character(Sys.time())
cat("\nStarting Forecast MCMC Sampler", tm, "\n" )


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
ss.mcmc = as.integer( args[2] )
thin.mcmc = as.integer( args[3] )
save( list=c("ss.mcmc", "thin.mcmc"), file=paste0(dir, "/src/mcmc_params.RData") )

## --------------------------------------------------------------------------------------
tm = as.character(Sys.time())
cat("\nLoading data", tm , "\n" )

unlargest200 = fread(file=paste0(dir, "/data/200isoRegionCodes.csv"))[order(iso)]
isos = unlargest200[order(iso)][,iso]

abelCohenTrain =
  fread(
    file=paste0(dir, "/data/abelCohen2019flowsv6_flowdt.csv")
    )[orig %in% isos & dest %in% isos][,`:=`(ac2019=flow, flow=NULL)][]


azoseRaftery2019 =
  fread(
    file=paste0(dir, "/data/azoseRaftery2019flows.csv" )
    )[,.(orig_code=origin, 
         dest_code=destination, 
		 orig=origIso,
         dest=destIso,
         year0=year,
         flow=migrantCount)][
           orig %in% isos & dest %in% isos][,`:=`(ar2019=flow, flow=NULL)][]


flowTrain =
  merge( abelCohenTrain, azoseRaftery2019, 
         by=c("year0", "orig", "dest"), all=TRUE)[
           ,.(flow=ifelse(is.na(ac2019), ar2019, ac2019)), 
           by=.(year0, orig, dest)] %>%
  merge( unique( azoseRaftery2019[,.(orig, dest, orig_code, dest_code)] ), 
         by=c("orig", "dest"))


# marginal flows
marginflowTrain = 
  merge(flowTrain[,.(inflow=sum(flow)), .(year0, iso=dest)], 
        flowTrain[,.(outflow=sum(flow)), .(year0, iso=orig)], 
        by=c("year0", "iso"))[,netflow:=inflow-outflow]

periodTrain = flowTrain[, sort( unique( year0 ) ) ]

nt = length( periodTrain )
nc = length( isos )


## --------------------------------------------------------------------------------------

# Population
data("pop", package="wpp2019")
data("popproj", package="wpp2019")

wppHistoricPop = as.data.table(pop)
wppProjectedPop = as.data.table(popproj)

popdt =
  merge(wppHistoricPop, wppProjectedPop, by=c("name", "country_code")) %>%
  merge(unlargest200, by="country_code")  %>%
  setnames(old="name", new="country_name") %>%
  setcolorder(c("area_code", "reg_code", "country_code",
                "area_name", "reg_name", "country_name", "iso"))

popNames = c("iso", as.character(seq(1950, 2020, 5)))
popIsodt =
  data.table::melt(popdt[,eval(popNames), with=FALSE],
       id.vars="iso", variable.name="year0", value.name="midPeriodPop000",
       variable.factor = FALSE)[,year0:=as.numeric(year0)][iso %in% isos][
         ,`:=`(pop=midPeriodPop000*1e3, midPeriodPop000=NULL)
       ][]


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


tmp = 
  marginflowTrain[,.(iso, year0, outflow)] %>% 
  dcast(formula=iso~year0, value.var="outflow")

outflowTrain = 
  tmp[,-"iso"] %>%
  as.matrix() %>%
  set_rownames(tmp$iso)

############################### BHM FORECAST ######################################

# init values
muHat = rowMeans( log( outrateTrain ) )

mu0 = mean( muHat )
sigma0 = sd( muHat )


a0 = 5.318086
b0 = 9.068054

# model
outflowCode <- nimbleCode({

  #### DEPART RATE MODEL ####
  for(c in 1:nc){
    for(t in 1:(nt-1)){
      y[c,t] ~ dnorm( mean=rho[c,t], sd=sigma[c] )
      rho[c,t] <- (1-phi) *mu[c] + phi * ylast[c,t]
    }
    
    ### GLOBAL PARAMETERS ###
    mu[c] ~ dnorm( mean=tau, sd=sigma0/sqrt(k0) ) 
    sigma[c] ~ dbeta(a0, b0)
  }
  
  phi ~ dunif(0,1)
  
  ### HYPER PRIORS ###
  tau ~ dnorm( mean=mu0, sd=100 )
})

outflow <- nimbleModel(code = outflowCode,
                    name = "outflow",
                    constants = list(mu0=mu0,
                                     sigma0=sigma0, 
                                     k0=1/10,
                                     a0=a0,
                                     b0=b0
                                     ),
                    data = list( y=log( outrateTrain )[,2:nt], 
                                 ylast=log( outrateTrain )[,1:(nt-1)]),
                    inits = list(mu=muHat, 
                                 sigma=rbeta(nc,a0,b0), 
                                 tau=mu0,
                                 phi=0.5
                                 ),
                    calculate=FALSE)

outflowConf <- configureMCMC(model=outflow, 
                          monitors=list("mu", "sigma", "phi"),
                          thin=10,
                          useConjugacy=FALSE,
                          print=FALSE)

Coutflow <- compileNimble(outflow)
outflowMCMC <- buildMCMC(outflowConf)
CoutflowMCMC <- compileNimble(outflowMCMC, project=outflow, resetFunctions=TRUE)

################################# BHM SAMPLE #########################################

tm = as.character(Sys.time())
cat("\nStarting sampler at", tm, "\n" )

CoutflowMCMC$run( ss.mcmc )
postSample = CoutflowMCMC$mvSamples %>% as.matrix() %>% as.data.table()
thin.ss = nrow( postSample )

origCodes = unlargest200[order(iso)][,.(i=1:.N, country_code)]

# mu
mu.mcmc = 
  postSample[,.SD, .SDcols=which( str_detect(colnames(postSample), "mu\\["))] %>%
  t() %>%
  as.data.table(keep.rownames=TRUE) %>%
  set_colnames(c("i", 1:thin.ss)) %>%
  .[,i:=tstrsplit(i, "\\[", keep=2)] %>%
  .[,i:=tstrsplit(i, "\\]", keep=1)] %>%
  .[,i:=as.numeric(i)] %>%
  .[order(i)] %>%
  merge(origCodes, by="i") %>%
  setcolorder(c("country_code", "i")) %>%
  .[order(i)] %>%
  .[,i:=NULL] %>%
  .[]
fwrite(mu.mcmc, file=paste0(dir, "/forecast/mu.csv"))


# sigma
sigma.mcmc =
  postSample[,.SD, .SDcols=which( str_detect(colnames(postSample), "sigma\\["))] %>%
  t() %>%
  as.data.table(keep.rownames=TRUE) %>%
  set_colnames(c("i", 1:thin.ss)) %>%
  .[,i:=tstrsplit(i, "\\[", keep=2)] %>%
  .[,i:=tstrsplit(i, "\\]", keep=1)] %>%
  .[,i:=as.numeric(i)] %>%
  .[order(i)] %>%
  merge(origCodes, by="i") %>%
  setcolorder(c("country_code", "i"))%>%
  .[order(i)] %>%
  .[,i:=NULL] %>%
  .[]
fwrite(sigma.mcmc, file=paste0(dir, "/forecast/sigma.csv"))


# phi
phi.mcmc =  postSample[, .(phi)]
fwrite(phi.mcmc, file=paste0(dir, "/forecast/phi.csv"))

tm = as.character(Sys.time())
cat("\nOutflow sampler completed at", tm, "\n" )

