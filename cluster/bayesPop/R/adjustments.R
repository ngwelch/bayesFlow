adjust.trajectories <- function(country, env, quant.env, adj.env=NULL) {
	if(is.null(adj.env)) adj.env <- new.env()
	wpp.year <- quant.env$wpp.year
	datasets <- list(totp='', totp.hch='', totpf='Fage', totpm='Mage', totpm.hch='Mage', totpf.hch='Fage')
	country.char <- as.character(country)
	for(traj.name in names(datasets)) {
		adj.name <- datasets[[traj.name]]
		dif.name <- paste0('AdjDpop', adj.name)
		
		if(is.null(adj.env[[dif.name]])) {
			#print(c(dif.name, adj.name))
			q <- quant.env[[paste0('quantiles', adj.name)]]
			adjust.quantiles(q, adj.name, wpp.year=wpp.year, env=adj.env)
		}
		dif <- if(length(dim(adj.env[[dif.name]]))>2) adj.env[[dif.name]][country.char,,] else adj.env[[dif.name]][country.char,,drop=FALSE]
		res <- env[[traj.name]]		
		if(length(dim(res))>2) { # includes age
			res21 <- aaply(res[1:21,,], 3, '-', dif)
			res21 <- aperm(res21, c(2,3,1))
			res[1:21,,] <- res21
		} else {
			res <- aaply(res, 2, '-', dif)
			res <- aperm(res, c(2,1))
		}
		if(is.list(res)) stop('')
		#if(traj.name == 'totpf') stop('')
		env[[traj.name]] <- res
	}
}

adjust.quantiles <- function(q, what, wpp.year, env=NULL) {
	dif <- NULL
	if(!is.null(env)) {
		if(!is.null(env[[paste0('AdjQpop', what)]])) return(env[[paste0('AdjQpop', what)]])
		if(!is.null(env[[paste0('AdjDpop', what)]])) dif <- env[[paste0('AdjDpop', what)]]
	}
	if(is.null(dif)) {
		if(is.null(env)) env <- new.env()
		countries <- dimnames(q)[[1]]
		ages <- NULL
		if(length(dim(q))>3) { # includes age dimension
			ages <- dimnames(q)[[2]]
			ages <- ages[as.numeric(ages)<=100]
		}
		wpp <- .get.wpp(env, what, countries, ages, wpp.year=wpp.year)
		if(length(dim(q))>3) { # includes age dimension
			years <- as.numeric(dimnames(q)[[4]])
			if((years[1] %% 5) != 0) years <- years+2 
			med.raw <- q[,,'0.5',as.character(years)%in%dimnames(wpp)[[3]]]
			if(length(dim(med.raw))>2) { # multiple countries
				med <- med.raw[,1:21,] # collapse to 21 age categories
				med[,21,] <- med.raw[,21,] + apply(med.raw[,22:27,], c(1,3), sum) 
			} else { #1 country
				med <- med.raw[1:21,]
				med[21,] <- med.raw[21,] + apply(med.raw[22:27,], 2, sum) 
				med <- abind(med, along=0) # add dimension
			}
			dif <- abind(matrix(0, nrow=dim(med)[1], ncol=21), med-wpp, along=3)		
		} else {
			years <- as.numeric(dimnames(q)[[3]])
			if((years[1] %% 5) != 0) years <- years+2 
			med <- q[,'0.5',as.character(years)%in%colnames(wpp)]
			dif <- as.matrix(cbind(0, med-wpp))
		}
	} else countries <- dimnames(dif)
	if(length(dim(q))>3) {
		if(dim(q)[1]==1){ # one country - the generic aaply fails because of some dimension dropping 
			res21 <- aaply(q[,1:21,,], 2, '-', dif[1,,], .drop=FALSE)
			res21 <- aperm(res21, c(2,1,3))
		} else {
			res21 <- aaply(q[,1:21,,], 3, '-', dif, .drop=FALSE)
			res21 <- aperm(res21, c(2,3,1,4))
		}
		res <- q
		res[,1:21,,] <- res21
	} else { # no age dimension
		res <- aaply(q, 2, '-', dif, .drop=FALSE)
		res <- aperm(res, c(2,1,3))
	}
	if(is.null(dimnames(dif)[[1]])) dimnames(dif)[[1]] <- countries
	if(is.null(dimnames(res)[[1]])) dimnames(res)[[1]] <- countries
	env[[paste0('AdjDpop', what)]] <- dif
	env[[paste0('AdjQpop', what)]] <- res
	return(res)
}

.get.wpp <- function(env, what, countries=NULL, ages=NULL, ...) {
	switch(which(c('', 'M', 'F', 'Mage', 'Fage') == what), 
				tpop(countries, prediction.only=TRUE, e=env, ...),
				tpopM(countries, prediction.only=TRUE, e=env, ...),
				tpopF(countries, prediction.only=TRUE, e=env, ...),
				tpopM(countries, prediction.only=TRUE, sum.over.ages=FALSE, ages=ages, e=env, ...),
				tpopF(countries, prediction.only=TRUE, sum.over.ages=FALSE, ages=ages, e=env, ...)
			)
}

if.not.exists.load <- function(name, env, wpp.year=2012) {
	if(!exists(name, where=env, inherits=FALSE)) {
		do.call('data', list(name, package=paste0('wpp', wpp.year), envir=env))
	    env[[name]] <- as.data.table(env[[name]])
	}
}

tpop <- function(countries, prediction.only=FALSE, e=NULL, ...) {
	# Create a dataset of total population
	if(is.null(e)) e <- new.env()
	if(!prediction.only) {
		if.not.exists.load('popM', e, ...)
		if.not.exists.load('popF', e, ...)
		tpop.obs <- sumMFbycountry('popM', 'popF', e)
	}
	#projection stored separately from observations
	if.not.exists.load('popMprojMed', e, ...)
	if.not.exists.load('popFprojMed', e, ...)
	tpopp <- sumMFbycountry('popMprojMed', 'popFprojMed', e)
	if(!prediction.only) tpopp <- merge(tpop.obs, tpopp, by='country_code')
	return(.reduce.to.countries(tpopp, countries))
}

tpopF <- function(...) return(tpop.sex('F', ...))
tpopM <- function(...) return(tpop.sex('M', ...))

tpop.sex <- function(sex, countries, sum.over.ages=TRUE, ages=NULL, prediction.only=FALSE, e=NULL, ...) {
	# Create a dataset of total population by sex
	if(is.null(e)) e <- new.env()
	if(!prediction.only) {
		dataset <- paste0('pop', sex)
		if.not.exists.load(dataset, e, ...)
		#do.call('data', list(dataset, package='wpp2012', envir=e))
		pop.obs <- if(sum.over.ages) sum.by.country(dataset) else sum.by.country.and.age(dataset)
	}
	dataset <- paste0('pop', sex, 'projMed')
	if.not.exists.load(dataset, e, ...)
	popp <- if(sum.over.ages) sum.by.country(e[[dataset]]) else sum.by.country.and.age(e[[dataset]])
	if(!prediction.only)  popp <- merge(pop.obs, popp, by='country_code')
	if(sum.over.ages) return(.reduce.to.countries(popp, countries))
	.reduce.to.countries.and.ages(popp, countries, ages)
}

.reduce.to.countries <- function(dataset, countries){
    tpop <- as.data.frame(dataset)
	tpop <- tpop[,-which(colnames(dataset)=='country_code')]
	rownames(tpop) <- dataset$country_code
	tpop[countries,]
}

.reduce.to.countries.and.ages <- function(dataset, countries, ages){
	dataset <- as.data.frame(dataset[dataset$country_code %in% as.integer(countries),])
	if(is.null(ages)) ages <- as.character(seq(0,100, by=5))
	age.vector <- as.character(dataset$age[1:21])
	age.vector <- unlist(strsplit(gsub('\\+', '-130', age.vector), '-'))
	age.vector <- age.vector[seq(1,length(age.vector), by=2)]
	colidx <- (1:ncol(dataset))[-which(colnames(dataset) %in% c('country_code', 'age'))]
	res <- array(NA, c(length(countries), length(ages), ncol(dataset)-2))
	for(i in 1:length(countries)) {
		idx <- which(dataset$country_code==countries[i])
		if(length(idx) == 0) next
		tmp <- dataset[idx,colidx]
		rownames(tmp) <- age.vector
		res[i,,] <- as.matrix(tmp[ages,])
	}
	dimnames(res)<- list(countries, ages, colnames(dataset)[colidx])
	res
}

sum.by.country <- function(dataset) {
	year.cols <- grep('^[0-9]{4}', colnames(dataset), value = TRUE)
	dataset[, c("country_code", year.cols), with = FALSE][, lapply(.SD, sum, na.rm = TRUE), by = "country_code"]
}

sum.by.country.and.age <- function(dataset) {
	year.cols <- grep('^[0-9]{4}', colnames(dataset), value = TRUE)
	dataset[, c("country_code", "age", year.cols), with = FALSE][, lapply(.SD, sum, na.rm = TRUE), by = c("country_code", "age")]
}

sumMFbycountry <- function(datasetM, datasetF, e) {
	tpopM <- sum.by.country(e[[datasetM]])
	tpopF <- sum.by.country(e[[datasetF]])
	tpopM[, 2:ncol(tpopM)] <- tpopM[,2:ncol(tpopM)] + tpopF[,2:ncol(tpopF)]
    tpopM
}

adjust.to.dataset <- function(country, q, adj.dataset=NULL, adj.file=NULL, years=NULL, use=c('write', 'trajectories')) {
	if(is.null(adj.dataset)) {
		adj.dataset <- read.table(adj.file, header=TRUE, check.names=FALSE)
	}
	colidx <- if(is.null(years)) (1:ncol(adj.dataset))[-which(colnames(adj.dataset)%in%c('country_code', 'country', 'name'))] else as.character(years)
	idx1 <- which(adj.dataset$country_code == country)
	if(use=='write') {
		med <- q['0.5']
		dif <- med - adj.dataset[idx1,colidx]
		return(q-dif)
	}
	if(use=='trajectories') {
		med <- apply(q, 1, 'median')[colnames(adj.dataset[,colidx])]
		dif <- as.matrix(med - adj.dataset[idx1,colidx])
		res <- aaply(q[colnames(adj.dataset[,colidx]),], 2, '-', dif)
		res <- aperm(res, c(2,1))
		if(! rownames(q)[1] %in% rownames(res))
			res <- rbind(q[1,], res) # add current year
		rownames(res) <- rownames(q)
		return(res)
	}
	return(NULL)
}

adjust.migration.if.needed <- function(time, year, country.codes, inputs, env) {
    # Adjust migration trajectories to given datasets
    country.codes.char <- as.character(country.codes)
    adj.names <- list(Male = "migration.adjustM.to", Female = "migration.adjustF.to", Total = "migration.adjust.to")
    period <- paste(year, year + 5, sep = "-")
    for(icntry in 1:length(country.codes)) {
        inp <- inputs[[country.codes.char[icntry]]]
        adjtype <- intersect(unlist(adj.names), names(inp))
        if(length(adjtype) == 0) next
        adjust <- list(Male = adj.names$Male %in% adjtype, Female = adj.names$Female %in% adjtype)
        adjust$Total <- !adjust$Male && !adjust$Female && adj.names$Total %in% adjtype
        cntry <- country.codes[icntry]
        # current migration projections
        migMage <- env$migm[,country.codes.char[icntry],]
        migFage <- env$migf[,country.codes.char[icntry],]
        migM <- colSums(migMage)
        migF <- colSums(migFage)
        lage <- nrow(migMage)
        current.mig.age <- list(Male = migMage, Female = migFage, Total = migMage + migFage)
        current.mig <- list(Male = migM, Female = migF, Total = migM + migF)
        
        # medians of the projections
        mean.mig <- age.props <- NULL
        for(whatadj in names(current.mig)) 
            mean.mig[[whatadj]] <- mean(current.mig[[whatadj]])
 
        adjusted <- NULL
        for(whatadj in names(adjust)) {
            if(!adjust[[whatadj]]) next
            # determine adjustment value
            adjds <- inp[[adj.names[[whatadj]]]]
            adjust.value <- adjds[, period]
            dif <- adjust.value - mean.mig[[whatadj]]
            # distribute dif among ages
            if(whatadj == "Total") { # if total, distribute dif also among sexes
                for(sex in c("Male", "Female")) {
                    age.props <- current.mig.age[[sex]]/t(current.mig[["Total"]])[rep(1,lage),]
                    adjusted[[sex]] <- current.mig.age[[sex]] + age.props * dif 
                }
            } else {
                age.props <- current.mig.age[[whatadj]]/t(current.mig[[whatadj]])[rep(1,lage),]
                adjusted[[whatadj]] <- current.mig.age[[whatadj]] + age.props * dif 
            }
        }
        env$migm[,country.codes.char[icntry],] <- adjusted$Male
        env$migf[,country.codes.char[icntry],] <- adjusted$Female
    }
}
