hapSim <- function(nall, n.strata = 1, cases, controls, haplo.freq, RR, RRcm, RRcf, RRstar, RR.mat, RRstar.mat, gen.missing.cases = NULL, gen.missing.controls = NULL, n.sim = 1000, xchrom = F, sim.comb.sex = "double", BR.girls, dire = "simfiles", ask = TRUE, verbose = TRUE, cpus = 1){
##
##
##
## Simulates genetic data in Haplin format, consisting of fetal effects, maternal effects and/or parent-of-origin effects. 
## Allows for simulation of both autosomal and X-linked markers, assuming HWE.
## Enables simulation of environmental data, i.e the input (relative risks, number of cases etc) may vary across different strata.
##
## nall is a vector of the number of alleles at each locus.
## n.strata is the number of strata.
## cases is a list of the number of case families. Each element contains the number of families of a specified family design. The names of the elements, i.e. the family designs, must be "mfc" (full triad), "mc" (mother-child dyad), "fc" (father-child dyad) or "c" (a single case child). 
## controls is a list of the number of control families. Each element contains the number of families of a specified family design. The possible family designs are "mfc" (full triad), "mc" (mother-child-dyad), "fc" (father-child dyad), "mf" (mother-father dyad), "c" (a single control child), "m" (a single control mother) or "f" (a single control father), i.e., the name of each element must equal one of the options. 
## haplo.freq is a numeric vector of the haplotype frequencies.
## RR is a numeric vector of relative risks.
## RRcm and RRcf are numeric vectors of the relative risks associated with the haplotypes transmitted from the mother and father, respectively.
## RRstar is a numeric vector and estimates how much double-dose children would deviate from the risk expected in a multiplicative dose-response relationship.
## RR.mat and RRstar.mat have an interpretation simular to RR and RRstar respectively, when simulating genetic data with maternal effects.
## gen.missing.cases generates missing values at random for the case families. Set to "NULL" by default, i.e. no missing values generated.
## gen.missing.controls generates missing values at random for the control families. Set to "NULL" by default, i.e. no missing values generated.
## n.sim is the number of simulations, i.e., the number of simulated data files.
## xchrom. Equals "FALSE" by default, which indicates simulation of autosomal markers. If "TRUE", hapSim simulates X-linked genes.
## sim.comb.sex. To be used with xchrom = TRUE. A character value that specifies how to handle gender differences on the X-chromosome. If "single", the effect of a (single) allele in males is equal to the effect of a single allele dose in females, and similarly, if "double", a single allele in males has the same effect as a double allele dose in females. Default is "double", which corresponds to X-inactivation.
## BR.girls. To be used with xchrom = TRUE. Gives the ratio of baseline risk for females to the baseline risk for males.
## dire gives the directory of the simulated data files.
## ask is a logical variable. If "TRUE", hapSim will ask before overwriting the files in an already existing directory.
## verbose. Logical. Default is "TRUE", which means that the file name is displayed for each iteration.
## cpus. Allows parallel processing of its analyses. The cpus argument should preferably be set to the number of available cpu's. If set lower, it will save some capacity for other processes to run. Setting it too high should not cause any serious problems.
##
##
##
##
#
.sim.maternal <- FALSE
if(!missing(RR.mat) | !missing(RRstar.mat)) .sim.maternal <- TRUE
#
.sim.poo <- FALSE
if(!missing(RRcm) | !missing(RRcf)) .sim.poo <- TRUE
#
if(.sim.poo && !missing(RR)) stop("RR cannot be present at the same time as RRcm and RRcf", call. = F)
#
.n.strata <- n.strata
#
## The number of possible haplotypes
.nhaplo <- prod(nall)
#
## The number of loci
.nloci <- length(nall)
#
.missing.controls <- FALSE
if(missing(controls)) .missing.controls <- TRUE
if(.missing.controls) .controls <- c(mfc=0)
else .controls <- controls
#
if(xchrom & sim.comb.sex %in% c("females","males")) BR.girls <- 1
#
## Defining arguments for each stratum
.arg <- list(nall=nall, n.strata=n.strata, cases=cases, controls=.controls, haplo.freq=haplo.freq, gen.missing.cases=gen.missing.cases, gen.missing.controls=gen.missing.controls, sim.maternal=.sim.maternal, sim.poo=.sim.poo, xchrom=xchrom, sim.comb.sex=sim.comb.sex, nhaplo=.nhaplo, nloci=.nloci, n.sim=n.sim)
if(.sim.poo) .arg <- c(.arg, list(RRcm=RRcm, RRcf=RRcf, RRstar=RRstar))
else .arg <- c(.arg, list(RR=RR, RRstar=RRstar))
if(.sim.maternal) .arg <- c(.arg, list(RR.mat=RR.mat, RRstar.mat=RRstar.mat))
if(xchrom) .arg <- c(.arg, list(BR.girls=BR.girls))
.strat.arg <- do.call(f.hapArg,args=.arg)
#
## Misc tests
lapply(1:.n.strata, function(x){do.call(f.hapTests, args = .strat.arg[x, ])})
#
.f.hapSim.arg <- .strat.arg[, -which(colnames(.strat.arg)%in%c("cases","controls","n.sim")), drop=FALSE]
#
if(nzchar(dire)){
	## If directory exists & ask == TRUE, query user
	if(file.exists(dire) & length(list.files(dire)) != 0 & ask){
		.answer <- readline(paste('Do you really want to overwrite ALL files in directory ', dire, '? (y/n)', sep = ""))
		if(.answer != "y"){
			cat("Stopped without overwriting file(s) in directory\n")
			return(invisible(FALSE))
		}
	}
	#
	unlink(dire, recursive = TRUE)
	dir.create(dire)
}	
#
if(xchrom & sim.comb.sex == "males") cat("The males are simulated assuming no contribution from fathers to sons\n")
#
## Runs in parallel only if cpus is numeric > 1
if(!missing(cpus)){
	if(!is.numeric(cpus)) stop('The number of cpu-s "cpus" must be numeric!', call. = F)
	.run.para <- (cpus > 1)
} else .run.para <- FALSE
#
.filename <- paste(dire, "/sim", 1:n.sim, ".dat", sep = "")
.filename.R <- paste("sim", 1:n.sim, sep = "")
#	
## If pure case/contol design, delete columns corresponding to mothers and fathers
if(all(sapply(.strat.arg[,"case.design"], length)==1) & all(sapply(.strat.arg[,"control.design"], length)==1) && (all(sapply(.strat.arg[,"case.design"], function(x){x=="c"}) & all(sapply(.strat.arg[,"control.design"], function(x){x=="c"}))))){
	cat("The files contain only the case-control data; the columns related to the parents will be deleted\n")
	.colnames <- paste("l_", rep(1:.nloci,each=2), ".c", c(1,2), sep = "")
} else .colnames <- paste("l_", rep(1:.nloci,each=6), rep(c(".m",".f",".c"),each=2), c(1,2), sep = "")
#
if(xchrom) .colnames <- c("sex",.colnames)
if(!.missing.controls) .colnames <- c("cc",.colnames)
if(.n.strata > 1) .colnames <- c("strat", .colnames)
#
.simulations <- function(i){
	.return <- NA
	#
	if(verbose) cat("Simulating ",.filename.R[i], ": \n", sep = "")
	#
	## Haplin must be made available in each run
	if(verbose) on.exit(if(!is.data.frame(.return) && is.na(.return)){cat("Run failed\n")}, add = T)	
	#
	## Haplin must be made available in each run
	suppressPackageStartupMessages(loadNamespace("Haplin"))
	#
	suppressMessages(.hapSim <- lapply(1:.n.strata, function(x){
		.hapSim <- do.call(f.hapSim, args = .f.hapSim.arg[x,])
		.hapSim <- cbind(strata=x, .hapSim)}
	))
	#
	.file <- do.call("rbind", .hapSim)
	#
	## If pure case/contol design, delete columns corresponding to mothers and fathers
	## .l is the number of columns to the left of the genetic data
	.l <- 2
	if(xchrom) .l <- 3
	.cc.design <- FALSE
	if(all(sapply(.strat.arg[,"case.design"], length)==1) & all(sapply(.strat.arg[,"control.design"], length)==1) && (all(sapply(.strat.arg[,"case.design"], function(i){i=="c"}) & all(sapply(.strat.arg[,"control.design"], function(i){i=="c"}))))) .cc.design <- TRUE
	.mf.col <- rbind(seq(1+.l,6*.nloci+.l,6))
	.mf.col <- c(.mf.col, .mf.col+2)
	.mf.col <- sort(c(.mf.col,.mf.col+1))
	if(.cc.design).file <- .file[, -.mf.col]
	if(.missing.controls) .file <- .file[,-2]
	if(.n.strata == 1) .file <- .file[,-1]
	#
	names(.file) <- .colnames
	#
	## Write to file
	if(nzchar(dire)) write.table(.file, file = .filename[i], quote = F, col.names = F, row.names = F)
	#
	if(verbose) cat("OK\n")
	if(dire=="") .return <- .file
	else .return <- TRUE
	return(invisible(.return))
}
#
## Parallel
cat("\n--- Running hapSim using ", cpus, " cpu(s) ---\n\n", sep = "")
#
if(!.run.para){
	.sim <- lapply(seq(length.out = n.sim), function(x) try(.simulations(x), silent = T))
	.errs <- sapply(.sim, class) == "try-error"
} else{
	w <- makeCluster(spec = cpus)
	on.exit(stopCluster(w), add = T)
	.sim <- parLapply(cl = w, 1:n.sim, fun = function(x) try(.simulations(x), silent = T))
	.errs <- parSapply(w, .sim, function(x){identical(class(x), "try-error")})
}
#
names(.sim) <- .filename.R
#
## Go through possible errors
if(any(.errs)){
	## Collect error messages, change output to NA with err. mess. as attributes
	.mess <- .sim[.errs]
	.err.res <- lapply(.mess, function(x) {
		.tmpres <- NA
		attributes(x) <- NULL
		attr(.tmpres, "Error") <- x
		return(.tmpres)
	})
	#
	.sim[.errs] <- .err.res
	#
	cat("\n --- hapSim has completed with errors ---\n\n")
	#
	## Display error messages
	cat("Error message(s):\n")
	#
	.sim.error <- .sim[which(is.na(.sim))]
	for(i in 1:length(.sim.error)){
		cat(paste("Did not generate ", names(.sim.error[i]), ": ", attr(.sim.error[[i]], "Error"), sep = ""), "\n")
		if(nzchar(dire)) cat(paste("Summary:\n",n.sim-length(.sim.error)," of ", n.sim, " files have been written to directory \"", dire, "\".", sep = ""), "\n")
	}
} else{
	cat("\n --- hapSim has completed without errors ---\n\n")
	if(nzchar(dire)) cat(paste("Summary:\nAll files have been written to directory \"", dire, "\".", sep = ""), "\n")
}
#
if(nzchar(dire)) cat(paste("The order of the columns in the files is: \n", paste(.colnames[1:min(15,length(.colnames))],sep="", collapse = " "), if(length(.colnames)>15) " etc", sep = "" ), "\n")
#
## Return
if(!nzchar(dire)) return(.sim)
else return(invisible(.sim))
}
