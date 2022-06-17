#!/usr/bin/env Rscript

# Rscript --no-save generate_forecast.R 10 10 

args = commandArgs(TRUE)

library(data.table)
library(magrittr)
library(dplyr)
library(lubridate)
library(wpp2019)
library(bayesTFR)
library(bayesLife)
library(package="bayesPop", lib.loc=file.path( getwd(), "lib" ))

batchSize = as.integer( args[1] )
batch = as.integer( args[2] )

# batchSize * batch should not exceed 1000--the number of pi samples

compute.summary = FALSE
debug = TRUE

if(debug) options(error=quote(dump.frames("last.dump", TRUE)))

dumpFile <- "last.dump.rda"
if (file.exists(dumpFile)) {
  file.remove(dumpFile)
}

# BayesPop Prototype
workingDir = "/homes/nwelch/research/raftery/welchRaftery2022revision"
pppDir <- file.path( workingDir, "data/ppp")

tfrDir <- file.path( pppDir, "tfr/results20190427BigSmall")
e0Dir <- file.path( pppDir, "LE/resultsHIV201904251050adjBigSmall")
migDir <- file.path( pppDir, "Mig")


pop.male.file <- file.path( pppDir, 'inputs/PopMale.txt')
pop.female.file <- file.path( pppDir, 'inputs/PopFemale.txt')
pasfr.file <- file.path( pppDir, 'inputs/PASFR.txt')
sex.ratio.file <- file.path( pppDir, 'inputs/SexRatioBirth.txt')

migration.male.file <- file.path( pppDir, 'Pop/patrick20190517/inputs/migrationM.txt')
migration.female.file <- file.path( pppDir, 'Pop/patrick20190517/inputs/migrationF.txt')
mx.male.file <- file.path( pppDir, 'Pop/patrick20190517/inputs/MxM.txt')
mx.female.file <- file.path( pppDir, 'Pop/patrick20190517/inputs/MxF.txt')

# use vwBaseYear2019 updated param info
control.file <- file.path( pppDir, 'inputs/vwBaseYear2019_modified_2020_11_19.txt')
my.locations.file <- file.path( pppDir, 'Pop/patrick20190517/inputs/UNlocations.txt')


## Flow Calculations

# Population
data("pop", package="wpp2019")
data("popproj", package="wpp2019")

wppHistoricPop = as.data.table(pop)
wppProjectedPop = as.data.table(popproj)
unlargest200 = fread( file.path(workingDir, "data/200isoRegionCodes.csv") )[order(iso)]

popdt = 
  merge(wppHistoricPop, wppProjectedPop, by=c("name", "country_code")) %>%
  merge(unlargest200, by="country_code")  %>%
  setnames(old="name", new="country_name") %>%
  setcolorder(c("area_code", "reg_code", "country_code", 
                "area_name", "reg_name", "country_name", "iso"))

popNames = c("iso", as.character(seq(1950, 2040, 5)))
popIsodt = 
  data.table::melt(popdt[,eval(popNames), with=FALSE], 
                   id.vars="iso", variable.name="year0", value.name="midPeriodPop000", 
                   variable.factor = FALSE)[,year0:=as.numeric(year0)][
                     ,`:=`(pop=midPeriodPop000*1e3, midPeriodPop000=NULL)
                   ][]

flow2019 = fread(file=file.path( workingDir, "data/abelCohen2019flowsv6_flowdt.csv") )

flowObs = 
  merge( flow2019, unlargest200[,.(orig=iso, orig_code=country_code)], by="orig") %>%
  merge( unlargest200[,.(dest=iso, dest_code=country_code)], by="dest") %>%
  .[,.(orig_code, dest_code, pop.year=year0+5, flow)]

marginflow2019 = 
  merge(flow2019[,.(inflow=sum(flow)), .(year0, iso=dest)], 
        flow2019[,.(outflow=sum(flow)), .(year0, iso=orig)], 
        by=c("year0", "iso"))[,netflow:=inflow-outflow]


popflow2019= merge(marginflow2019, popIsodt)[
  ,`:=`(netrate=netflow/(5*pop/1e3), 
        outrate=outflow/(5*pop/1e3))]

outrate2019 = popflow2019[,.(year0, iso, outrate)]
outRateLast = outrate2019[year0==2015][order(iso)][,outrate]
names(outRateLast) = unlargest200[order(iso),country_code]

tmp = dcast(outrate2019, formula=iso~year0, value.var="outrate")[order(iso)]
outrateMatrix = 
  as.matrix( tmp[,-"iso"] ) %>%
  set_rownames(tmp$iso)
rm(tmp)

data("migration", package="wpp2019")

wppMigration = as.data.table(migration)

# net migration in thousands of people over 5 year period
netMigrantsdt = 
  merge(wppMigration, unlargest200, by=c("country_code")) %>%
  setnames(old="name", new="country_name") %>%
  setcolorder(c("area_code", "reg_code", "country_code", 
                "area_name", "reg_name", "country_name", "iso"))

migPeriodNames = c("iso", paste0( seq(1950, 2040, 5), "-", seq(1955, 2045, 5) ))
unNetMigdt = 
  data.table::melt(netMigrantsdt[,eval(migPeriodNames), with=FALSE], 
                   id.vars="iso", variable.name="period", value.name="periodMigration000")[
                     ,year0:=tstrsplit(period, "-")[[1]]
                   ][,year0:=as.numeric(year0)][
                     ,`:=`(netflow=1000*periodMigration000, periodMigration000=NULL)
                   ][]


isos = unlargest200$iso

mu.mcmc = fread(file=file.path( migDir, "mu.csv") )
sigma.mcmc = fread(file=file.path( migDir, "sigma.csv"))
phi.mcmc = fread(file=file.path( migDir, "phi.csv") )$phi


## BayesPop Run
pred <- pop.predict(
  countries=unlargest200[order(iso),country_code],
  output.dir = workingDir,
  nr.traj = 1000,
  batch=batch,
  batch.size=batchSize,
  end.year = 2045, start.year = 1950,
  wpp.year = 2019,
  present.year = 2020, #2020
  verbose = TRUE, 
  replace.output = TRUE,
  parallel=FALSE,
  lc.for.hiv = TRUE,
  lc.for.all = TRUE,
  inputs = 
    list(
      popM = pop.male.file       # age-specific male pop - must include column "2020"
      , popF = pop.female.file   # age-specific female pop - must include column "2020"
      , tfr.sim.dir = tfrDir
      #, tfr.file = "median_"
      , pasfr = pasfr.file
      , srb = sex.ratio.file
      , e0F.sim.dir = e0Dir
      #, e0F.file = "median_"
      , e0M.sim.dir = 'joint_'
      #, e0M.file = "median_"
      , mxM = mx.male.file
      , mxF = mx.female.file
      , migM = migration.male.file      # age-specific male migration if available
      , migF = migration.female.file    # age-specific female migration if available
      , patterns = control.file	       # control file
    ),
  fixed.mx = FALSE, 
  fixed.pasfr = FALSE,
  keep.vital.events = FALSE,
  use.migration.flow.model = TRUE, 
  compute.summary = compute.summary,
  migration.settings = list(mu=mu.mcmc, phi=phi.mcmc, sigma=sigma.mcmc, 
                            ini.out.rates=outRateLast, 
                            flowData=flowObs,
                            pi.dir=migDir, 
							bound.pop.drop=TRUE)
)
