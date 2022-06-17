#!/usr/bin/env Rscript

library(data.table)
library(magrittr)
library(dplyr)

workingDir = "/homes/nwelch/research/raftery/welchRaftery2022"

env=new.env()

flowQuant.time = marginQuant.time = globalMigrationPost.time = list()
countryAgePop.time = countryPop.time = worldPop.batch = list()
countryAgePopQuant.time = countryPopQuant.time = worldPopQuant.time = list()
batchIndex = 1:10
for(t in 1:5){
	cat("Processing time period", t, "\n")	
  flow.batch = countryAgePop.batch = countryPop.batch = worldPop.batch = list()
  for(b in batchIndex){
	cat("Processing batch", b, "at", format(Sys.time(), "%X"), "\n")	
    batch = paste0("batch_", b)
    flow.female.bt = paste0("flow_female_time_", t, ".rda")
    flow.male.bt = paste0("flow_male_time_", t, ".rda")
    pop.bt = paste0("pop_time_", t, ".rda")

    flow.female.file = file.path(workingDir, "bayesPop.output", batch, "raw", flow.female.bt )
    flow.male.file = file.path(workingDir, "bayesPop.output", batch, "raw", flow.male.bt )
    pop.file = file.path(workingDir, "bayesPop.output", batch, "raw", pop.bt)

    load(file.path( flow.female.file), env=env)
    load(file.path( flow.male.file), env=env)
    load(file.path( pop.file), env=env)

    env$flow.age = env$flow.age.female + env$flow.age.male
    env$flow.total =
      apply( env$flow.age, MARGIN=c("origin", "destination", "trajectory"), sum)


    flow.batch[[b]] =
      reshape2::melt( env$flow.total, value.name="flow" ) %>%
      as.data.table() %>%
      .[,trajectory:=trajectory + (b-1)*100] %>%
      .[origin!=destination]

    popf =
      reshape2::melt( env$totpf,
                      varnames=c("age", "country_code", "trajectory"),
                      value.name="popf") %>%
      as.data.table() %>%
	  .[,popf:=1e3*popf] %>%
      .[,trajectory:=trajectory + (b-1)*100] %>%
      .[]

    popm =
      reshape2::melt( env$totpm,
                      varnames=c("age", "country_code", "trajectory"),
                      value.name="popm") %>%
      as.data.table() %>%
	  .[,popm:=1e3*popm] %>%
      .[,trajectory:=trajectory + (b-1)*100] %>%
      .[]

    countryAgePop.batch[[b]] = merge(popf, popm)[,popt:=popf+popm][]
    countryPop.batch[[b]] =
      countryAgePop.batch[[b]][,.(popf=sum(popf), popm=sum(popm), popt=sum(popt)),
            .(country_code, trajectory)]
    worldPop.batch[[b]] =
      countryPop.batch[[b]][,.(popf=sum(popf), popm=sum(popm), popt=sum(popt)), trajectory]

  }

  # flows
  flowPost.time = data.table::rbindlist( flow.batch[batchIndex] )
  flowQuant.time[[t]] =
    flowPost.time[, {
      q = quantile(flow, c(0.05, 0.1, 0.5, 0.9, 0.95), na.rm=TRUE)
      list(
        year=seq(2023, 2043, 5)[t],
        q05 = q[1],
        q10 = q[2],
        q50 = q[3],
        q90 = q[4],
        q95 = q[5]
        )
        }, keyby=list(origin, destination)]

  # marginal flows
  outflowPost.time =
    flowPost.time[,.(outflow=sum(flow)),.(country_code=origin, trajectory)]
  inflowPost.time =
    flowPost.time[,.(inflow=sum(flow)),.(country_code=destination, trajectory)]
  marginflowPost.time =
    merge(outflowPost.time, inflowPost.time)[,netflow:=inflow-outflow][]


  marginQuant.time[[t]] =
    marginflowPost.time[, {
      qin = quantile(inflow, probs=c(0.05, 0.1, 0.5, 0.9, 0.95), na.rm=TRUE);
      qout = quantile(outflow, probs=c(0.05, 0.1, 0.5, 0.9, 0.95), na.rm=TRUE);
      qnet = quantile(netflow, probs=c(0.05, 0.1, 0.5, 0.9, 0.95), na.rm=TRUE);
      list(
        year=seq(2023, 2043, 5)[t],
        q05in = qin[1],
        q10in = qin[2],
        q50in = qin[3],
        q90in = qin[4],
        q95in = qin[5],
        q05out = qout[1],
        q10out = qout[2],
        q50out = qout[3],
        q90out = qout[4],
        q95out = qout[5],
        q05net = qnet[1],
        q10net = qnet[2],
        q50net = qnet[3],
        q90net = qnet[4],
        q95net = qnet[5]
        )
        }, country_code]

  # population
  countryAgePop.time = data.table::rbindlist( countryAgePop.batch[batchIndex] )
  countryAgePopQuant.time[[t]] =
    countryAgePop.time[, {
      qf = quantile(popf, probs=c(0.05, 0.1, 0.5, 0.9, 0.95), na.rm=TRUE);
      qm = quantile(popm, probs=c(0.05, 0.1, 0.5, 0.9, 0.95), na.rm=TRUE);
      qt = quantile(popt, probs=c(0.05, 0.1, 0.5, 0.9, 0.95), na.rm=TRUE);
      list(
        year=seq(2023, 2043, 5)[t],
        q05f = qf[1],
        q10f = qf[2],
        q50f = qf[3],
        q90f = qf[4],
        q95f = qf[5],
        q05m = qm[1],
        q10m = qm[2],
        q50m = qm[3],
        q90m = qm[4],
        q95m = qm[5],
        q05t = qt[1],
        q10t = qt[2],
        q50t = qt[3],
        q90t = qt[4],
        q95t = qt[5]
        )
        }, .(age, country_code) ]


  countryPop.time = data.table::rbindlist( countryPop.batch[batchIndex] )
  countryPopQuant.time[[t]] =
    countryPop.time[, {
      qf = quantile(popf, probs=c(0.05, 0.1, 0.5, 0.9, 0.95), na.rm=TRUE);
      qm = quantile(popm, probs=c(0.05, 0.1, 0.5, 0.9, 0.95), na.rm=TRUE);
      qt = quantile(popt, probs=c(0.05, 0.1, 0.5, 0.9, 0.95), na.rm=TRUE);
      list(
        year=seq(2023, 2043, 5)[t],
        q05f = qf[1],
        q10f = qf[2],
        q50f = qf[3],
        q90f = qf[4],
        q95f = qf[5],
        q05m = qm[1],
        q10m = qm[2],
        q50m = qm[3],
        q90m = qm[4],
        q95m = qm[5],
        q05t = qt[1],
        q10t = qt[2],
        q50t = qt[3],
        q90t = qt[4],
        q95t = qt[5]
        )
        }, country_code ]

  worldPop.time = data.table::rbindlist( worldPop.batch[batchIndex] )
  worldPopQuant.time[[t]] =
    worldPop.time[, {
      qf = quantile(popf, probs=c(0.05, 0.1, 0.5, 0.9, 0.95), na.rm=TRUE);
      qm = quantile(popm, probs=c(0.05, 0.1, 0.5, 0.9, 0.95), na.rm=TRUE);
      qt = quantile(popt, probs=c(0.05, 0.1, 0.5, 0.9, 0.95), na.rm=TRUE);
      list(
        year=seq(2023, 2043, 5)[t],
        q05f = qf[1],
        q10f = qf[2],
        q50f = qf[3],
        q90f = qf[4],
        q95f = qf[5],
        q05m = qm[1],
        q10m = qm[2],
        q50m = qm[3],
        q90m = qm[4],
        q95m = qm[5],
        q05t = qt[1],
        q10t = qt[2],
        q50t = qt[3],
        q90t = qt[4],
        q95t = qt[5]
        )
        }]
  
  # global flows
  globalMigrationPost.time[[t]] = 
		  merge( outflowPost.time[,.(flowt=sum(outflow)), trajectory],
				 worldPop.time[,.(popt, trajectory)], by="trajectory")[,{
					qpct=quantile( 100*flowt/popt, probs=c(0.05, 0.1, 0.5, 0.9, 0.95), na.rm=TRUE);
					qcnt=quantile( flowt, probs=c(0.05, 0.1, 0.5, 0.9, 0.95), na.rm=TRUE);
					list(
					year=seq(2023, 2043, 5)[t],
					q05pct = qpct[1],
					q10pct = qpct[2],
					q50pct = qpct[3],
					q90pct = qpct[4],
					q95pct = qpct[5],
					q05cnt = qcnt[1],
					q10cnt = qcnt[2],
					q50cnt = qcnt[3],
					q90cnt = qcnt[4],
					q95cnt = qcnt[5]
				)
		  }]
	
}

flowPost = data.table::rbindlist( flowQuant.time )
flowMarginPost = data.table::rbindlist( marginQuant.time )
globalMigrationPost = data.table::rbindlist( globalMigrationPost.time )

fwrite(flowPost, file=file.path(workingDir, "bayesPop.output", "flowPost.csv"))
fwrite(flowMarginPost, file=file.path(workingDir, "bayesPop.output", "flowMarginPost.csv"))
fwrite(globalMigrationPost, file=file.path(workingDir, "bayesPop.output", "globalMigrationPost.csv"))

countryAgePost =  data.table::rbindlist( countryAgePopQuant.time )
countryPost =  data.table::rbindlist( countryPopQuant.time )
worldPost = data.table::rbindlist( worldPopQuant.time )

fwrite(countryAgePost, file=file.path(workingDir, "bayesPop.output", "countryAgePost.csv"))
fwrite(countryPost, file=file.path(workingDir, "bayesPop.output", "countryPost.csv"))
fwrite(worldPost, file=file.path(workingDir, "bayesPop.output", "worldPost.csv"))



cat("Processing completed at", format(Sys.time(), "%X"), "\n")	
