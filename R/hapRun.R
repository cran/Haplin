hapRun <- function(nall, n.strata = 1, cases, controls, haplo.freq, RR, RRcm, RRcf, RRstar, RR.mat, RRstar.mat, hapfunc = "haplin", gen.missing.cases = NULL, gen.missing.controls = NULL, n.sim = 1000, xchrom = FALSE, sim.comb.sex = "double", BR.girls, dire, ask = TRUE, cpus = 1, slaveOutfile = "", ...){
##
##
## Simulates genetic data in Haplin format, consisting of fetal effects, maternal effects and/or parent-of-origin effects, then runs haplin or haplinSlides on the simulated data files. 
## Allows for simulations and calculations on both autosomal and X-linked markers, assuming HWE.
##
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
## hapfunc. Defines which haplin function to run, "haplin" or "haplinSlide,". "haplinStrat" is not yet implemented.
## gen.missing.cases generates missing values at random for the case families. Set to "NULL" by default, i.e. no missing values generated.
## gen.missing.controls generates missing values at random for the control families. Set to "NULL" by default, i.e. no missing values generated.
## n.sim is the number of simulations, i.e., the number of simulated data files.
## xchrom. Equals "FALSE" by default, which indicates simulation of autosomal markers. If "TRUE", hapSim simulates X-linked genes.
## sim.comb.sex. To be used with xchrom = TRUE. A character value that specifies how to handle gender differences on the X-chromosome. If "single", the effect of a (single) allele in males is equal to the effect of a single allele dose in females, and similarly, if "double", a single allele in males has the same effect as a double allele dose in females. Default is "double", which corresponds to X-inactivation.
## BR.girls. To be used with xchrom = TRUE. Gives the ratio of baseline risk for females to the baseline risk for males.
## dire gives the directory of the simulated data files. Missing by default, i.e. no files are written to directory
## ask is a logical variable. If "TRUE", hapSim will ask before overwriting the files in an already existing directory.
## cpus. Allows parallel processing of its analyses. The cpus argument should preferably be set to the number of available cpu's. If set lower, it will save some capacity for other processes to run. Setting it too high should not cause any serious problems. 
## slaveOutfile. Character. If slaveOutfile = "" (default), output from all running cores will be printed in the standard R session window.
## ... Arguments to be used by haplin or haplinSlides.
##
##
##
#
.sim.maternal <- FALSE
if(!missing(RR.mat) | !missing(RR.mat)) .sim.maternal <- TRUE
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
if(.missing.controls) controls <- c(mfc=0)
#
if(xchrom & sim.comb.sex %in% c("females","males")) BR.girls <- 1
#
## Defining arguments for each stratum
.arg <- list(nall=nall, n.strata=n.strata, cases=cases, controls=controls, haplo.freq=haplo.freq, gen.missing.cases=gen.missing.cases, gen.missing.controls=gen.missing.controls, sim.maternal=.sim.maternal, sim.poo=.sim.poo, xchrom=xchrom, sim.comb.sex=sim.comb.sex, nhaplo=.nhaplo, nloci=.nloci, n.sim=n.sim)
if(.sim.poo) .RR.arg <- list(RRcm=RRcm, RRcf=RRcf, RRstar=RRstar)
else .RR.arg <- list(RR=RR, RRstar=RRstar)
if(.sim.maternal) .RR.arg <- c(.RR.arg,list(RR.mat=RR.mat, RRstar.mat=RRstar.mat))
if(xchrom) .RR.arg <- c(.RR.arg,list(BR.girls=BR.girls))
.arg <- c(.arg, .RR.arg)
.strat.arg <- do.call(f.hapArg,args=.arg)
.f.hapSim.arg <- .strat.arg[, -which(colnames(.strat.arg)%in%c("cases","controls","n.sim")), drop=FALSE]
#
## Misc tests
lapply(1:.n.strata, function(x){do.call(f.hapTests, args = .strat.arg[x, ])})
#
## Misc errors
if(!missing(cpus) && !is.numeric(cpus)) stop('The number of cpu-s "cpus" must be numeric!', call. = F)
#
.hapfunc <- c("haplin", "haplinSlide", "haplinStrat")
if(!(hapfunc %in% .hapfunc)) stop("Argument \"hapfunc\" is not specified correctly", call. = F)
if(hapfunc == "haplin" & .n.strata != 1) stop("Function \"haplin\" is not valid when the number of strata is larger than 1", call. = F)
if(hapfunc == "haplinStrat" & .n.strata == 1) stop("Function \"haplinStrat\" is not valid when the number of strata is equal to 1", call. = F)
if(hapfunc == "haplinSlide" & .n.strata > 1) stop("The possibility to run \"haplinSlide\" when the number of strata is larger than 1 is not yet implemented", call. = F)
#
.dire <- FALSE
if(!missing(dire)) .dire <- TRUE
if(.dire){
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
## Choose Haplin design according to the arguments "cases" and "controls"
.design <- sapply(1:.n.strata, function(x){
	if(all(.strat.arg[,"control.mat"][[x]]==0) | .missing.controls){
		if(all(.strat.arg[,"case.design"][[x]]=="c")) stop("Only case children are given. No controls are available", call. = F)
		.design <- "triad"
	}	
	else if(all(.strat.arg[,"case.design"][[x]]=="c") & all(.strat.arg[,"control.design"][[x]]=="c")) .design <- "cc"
	else .design <- "cc.triad"
	return(.design)
})
if(length(unique(.design))==1) .design <- unique(.design)
else stop("Unable to specify haplin design due to the combination of arguments \"cases\" and \"controls\"", call.=F)
#
## Specify arguments to be passed to haplin, haplinStrat, haplinSlide
.n.vars <- 1
if(.missing.controls) .n.vars <- 0
.strata <- NULL
.ccvar <- NULL
.sex = NULL
#
if(.n.strata > 1){ 
	.n.vars <- .n.vars + 1
	.strata <- 1
}	
if(xchrom) .n.vars <- .n.vars + 1
if(.design != "triad" & !xchrom) .ccvar <- .n.vars
else if(.design == "triad" & xchrom) .sex <- .n.vars	
else if(.design != "triad" & xchrom){
	.ccvar <- .n.vars - 1
	.sex <- .n.vars
}
#
## SET UP TEMPORARY FILE FOR HAPLOTYPES
.tmphaplofile <- tempfile(tmpdir = ".")
on.exit(unlink(.tmphaplofile), add = T)
#
## WRITE (TEMPORARY) FILE CONTAINING HAPLOTYPES
.alleles <- lapply(nall,seq)
if(hapfunc=="haplinSlide"){
	if(any(nall!=2)) stop("\"haplinSlide\" is currently only implemented for diallelic loci", call. = F)
	.alleles <- lapply(rep(2,list(...)$winlength),seq)
	warning("\"haplinSlide\" is only partially implemented, and the results should be interpreted with caution. Please see the help file for more information", call. = F)
	.ref.cat <- 1
}	
.haplotypes <- do.call("expand.grid", .alleles)
.haplotypes <- f.create.tag(.haplotypes, sep = "-")
write.table(dframe(haplos = .haplotypes), file = .tmphaplofile, quote = F, row.names = F, col.names = T)
#
## FORCE haplin LATER ON TO USE SAME HAPLOTYPES AND SAME REFERENCE CATEGORY
.haplo.file <- .tmphaplofile
if(hapfunc != "haplinSlide") {
	if(is.list(haplo.freq)) .haplo.freq <- haplo.freq[[which.max(unlist(lapply(haplo.freq,max)))]]
	else .haplo.freq <- haplo.freq
	.ref.cat <- which.max(.haplo.freq)
}	
#
## Arguments specified by user to be passed on to by haplin, haplinSlide or haplinStrat
.lu <- list(...)
.ld <- list(n.vars = .n.vars, design = .design, ccvar = .ccvar, xchrom = xchrom, sex = .sex, strata = .strata, reference = .ref.cat, verbose = FALSE, haplo.file = .haplo.file, use.missing = TRUE, printout = FALSE)
.nu <- names(.lu)
.nd <- names(.ld) 
#
## Set argument(s) to user defined value(s)
.ld[.nu] <- .lu
#
## Report if default haplin arguments are overwritten
.fixed <- c("n.vars", "ccvar", "sex", "strata", "max.haplos", "haplo.file", "data.out")
if(any(.fixed%in%.nu)) stop(paste("The following arguments are not valid:", paste(intersect(.fixed,.nu), collapse = ", "), "\n", sep = " "), call. = F)
if(!is.numeric(.ld$reference) | .ld$reference>.nhaplo |.ld$reference<1) stop("Argument \"reference\" is not valid", call. = F)
.inter <- intersect(.nu,.nd)
if(length(.inter) != 0) invisible(cat(paste("The following arguments are overwritten by user:", paste(.inter, collapse = ", "), "\n", sep = " ")))
#
if(.dire) .simfiles <- paste(dire, "/", 1:n.sim, sep = "")
else .simfiles <- 1:n.sim
.filename <- paste(.simfiles, "/", "sim1.dat", sep = "")
#
### Function running haplin (one version or another)
.hapRun <- function(x){
	## Haplin must be made available in each run
	suppressPackageStartupMessages(require(Haplin, quietly = T))
	#
	if(.dire) cat("\nSimulations using ", .simfiles[x], "\n", sep = "")
	else cat("\nSimulations using file number ", .simfiles[x], "\n", sep = "")
	#
	suppressMessages(.hapSim <- lapply(1:.n.strata, function(i){
			.hapSim <- do.call(f.hapSim, args = .f.hapSim.arg[i,])
			.hapSim <- cbind(strata=i, .hapSim)}
	))
	#
	.file <- do.call("rbind", .hapSim)
	#
	## If pure case/contol design, delete columns corresponding to mothers and fathers
	## .l is the number of columns to the left of the genetic data
	.l <- 2
	if(xchrom) .l <- 3
	.mf.col <- rbind(seq(1+.l,6*.nloci+.l,6))
	.mf.col <- c(.mf.col, .mf.col+2)
	.mf.col <- sort(c(.mf.col,.mf.col+1))
	if(.design=="cc") .file <- .file[,-.mf.col]
	if(.missing.controls) .file <- .file[,-2]
	if(.n.strata == 1) .file <- .file[,-1]
	.file <- as.matrix(.file)
	#
	if(.n.vars > 0){
		if(.design=="cc") colnames(.file) <- c(colnames(.file)[1:.n.vars],paste("l_",rep(1:.nloci,each=2),rep(".c",each=2),c(1,2),sep=""))
		else colnames(.file) <- c(colnames(.file)[1:.n.vars],paste("l_",rep(1:.nloci,each=6),rep(c(".m",".f",".c"),each=2),c(1,2),sep=""))
	} else{
		if(.design=="cc") colnames(.file) <- paste("l_",rep(1:.nloci,each=2),rep(".c",each=2),c(1,2),sep="")
		else colnames(.file) <- paste("l_",rep(1:.nloci,each=6),rep(c(".m",".f",".c"),each=2),c(1,2),sep="")
	}	
	#
	.ld <- c(list(data=.file), allele.sep = " ", .ld)
	#
	## Write to file if .dire = T
	if(.dire && nzchar(dire)){
		dir.create(.simfiles[x])
		write.table(.file, file = .filename[x], quote = F, col.names = F, row.names = F)
	}
	#
	if(hapfunc == "haplin"){
		.haplin <- do.call(haplin, args = .ld)
		.output <- haptable(.haplin)
	}
	if(hapfunc == "haplinSlide"){
		.ld <- c(.ld, table.output = F)
		.haplin <- do.call(haplinSlide, args = .ld)
		.output <- suest(.haplin)
	}
	if(hapfunc == "haplinStrat"){
		.haplin <- do.call(haplinStrat, args = .ld)
		.output <- gxe(.haplin)
	}
	#
	return(invisible(.output))
}
#
## Parallel
w <- makeCluster(spec = cpus, outfile = slaveOutfile)
on.exit(stopCluster(w), add = T)
#
## Output
cat("\n--- Running hapRun using ", cpus, " cpu(s) ---\n", sep = "")
if(slaveOutfile != "") cat("\nOutput is written to \"", slaveOutfile, "\"\n", sep = "")
#
## Run .hapRun in parallel
.run <- parLapply(cl = w, 1:length(.simfiles), fun = function(x) try(.hapRun(x), silent = T))
names(.run) <- .simfiles
#
## Check for failures in function .hapRun
.errs <- parSapply(w, .run, function(x){identical(class(x), "try-error")})
if(any(.errs)){
	## Collect error messages, change output to NA with err. mess. as attributes
	.mess <- .run[.errs]
	.err.res <- lapply(.mess, function(x){
		.tmpres <- NA
		attributes(x) <- NULL
		attr(.tmpres, "Error") <- x
		return(.tmpres)
	})
	#
	.run[.errs] <- .err.res
	#
	cat(paste("\n--- hapRun has completed with errors ---\nhapRun was completed for ", length(.simfiles)-length(which(.errs==1)), " out of ", length(.simfiles), " simulated files \n\n", sep = ""))
	#
	## Display error messages
	cat("Error message(s): \n")
	.run.error <- .run[which(is.na(.run))]
	for(i in 1:length(.run.error)){
		cat(paste("hapRun failed for simulation ", names(.run.error[i]), ": ", attr(.run.error[[i]], "Error"), sep = ""))
	}
}else cat(paste("\nhapRun ran without errors. All ", n.sim, " simulations succeeded\n", sep = ""))	
#
## Stack dataframes from haplin runs into a single dataframe
if(hapfunc == "haplin"){
	.run <- toDataFrame(.run)
	colnames(.run)[1] <- "sim.no"
}
#
## Add hapfunc as attribute to .run
attr(.run,"hapfunc") <- hapfunc
#
## Return
return(.run)
}