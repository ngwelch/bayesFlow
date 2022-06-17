#!/usr/bin/env Rscript

library(data.table)
library(magrittr)
library(dplyr)
library(stringr)

args = commandArgs(TRUE) 
dir = as.character( args[1] )

## --------------------------------------------------------------------------------------

unlargest200 = fread("data/200isoRegionCodes.csv")[order(iso)]
countryCodes = unlargest200[order(iso)][,country_code] %>% as.integer()

countryCodes.kappa = 
		list.files( file.path( dir, "kappa") ) %>%
		strsplit(split=".", fixed=TRUE) %>%
		sapply( "[", 1 ) %>%
		as.integer()

countryCodes.psi = 
		list.files( file.path( dir, "psi") ) %>%
		strsplit(split=".", fixed=TRUE) %>%
		sapply( "[", 1 ) %>%
		as.integer()

countryCodes.success = intersect( countryCodes.kappa, countryCodes.psi )

# figure out which countries are missing
index.fail = which(!(countryCodes %in% countryCodes.success))

fileConn = file.path(dir, "sampler_fail.txt")
writeLines( as.character(index.fail), con=fileConn, sep="," )
close(fileConn)
