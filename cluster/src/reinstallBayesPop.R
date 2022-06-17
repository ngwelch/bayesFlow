#!/usr/bin/env Rscript

workingDir = "/homes/nwelch/research/raftery/welchRaftery2022revision"

if( ("bayesPop" %in% rownames(installed.packages())) ) remove.packages("bayesPop", lib=file.path( workingDir, "lib") )

install.packages(pkgs=file.path( workingDir, "bayesPop" ),
				 lib=file.path( workingDir, "lib" ) ,
				 repos=NULL,
				 dependencies=FALSE,
				 type="source",
				 clean=TRUE,
				 verbose=FALSE)


