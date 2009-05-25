haplin <- function(filename, 
markers = "ALL", n.vars = 0, sep = " ", allele.sep = ";", na.strings = "NA",
design = "triad", use.missing = FALSE, xchrom = FALSE, maternal = FALSE, scoretest = "no",
ccvar = NULL, covar = NULL, sex = NULL,
reference = "reciprocal", response = "free", threshold = 0.01, max.haplos = NULL, haplo.file = NULL,
resampling = FALSE, max.EM.iter = 50, data.out = FALSE, verbose = TRUE, printout = TRUE
)
#"haplin"<- function(filename, filespecs, model, variables, haplos, control)
{
##
## filename: (character string) filename 
## filespecs: (list) specifics of the file. n.vars is the number of columns before genetic data columns, sep = " " is the column separator, allele.sep = ";" separates alleles within marker, na.strings = "NA" identifies missing values
## variables: (list) ccvar is the position of the case-control variable, covar is the position of the environment covariate
## 
# PAABEGYNT 13/1-04 KL. 22.46
#
#
.mcall <- lapply(match.call()[-1], function(x) eval.parent(x, 3))
.defaults <- formals()
.info <- f.check.pars(.mcall, .defaults)

#
## ARGUMENTS NOT GIVEN ARE SET TO NULL (THEN LATER TO DEFAULT VALUES BY f.check.pars)
#if(missing(filename)) filename <- NULL
#if(missing(filespecs)) filespecs <- NULL
#if(missing(model)) model <- NULL
#if(missing(variables)) variables <- NULL
#if(missing(haplos)) haplos <- NULL
#if(missing(control)) control <- NULL
#
## COLLECT INFO, SET DEFAULT VALUES AND CHECK FOR SANITY
#.info <- list(filename = filename, filespecs = filespecs, model = model, variables = variables, haplos = haplos, control = control)
#.info <- f.check.pars(.info)
#
## SET PARAMETERS, FOR SIMPLICITY
design <- .info$model$design
xchrom <- .info$model$xchrom
maternal <- .info$model$maternal
max.EM.iter <- .info$control$max.EM.iter
#
# VET IKKE MED DISSE:
use.missing <- .info$model$use.missing
verbose <- .info$control$verbose
resampling <- .info$control$resampling
data.out <- .info$control$data.out
printout <- .info$control$printout

#
## MAKE SURE CONTRASTS ARE RIGHT (ALTHOUGH NOT STRICTLY NECESSARY):
.old.options <- options()
on.exit(options(.old.options))
options(contrasts = c("contr.treatment", "contr.poly"))	#
options(stringsAsFactors = F)
#
## INSTALL MASS (FOR THE mvrnorm FUNCTION):
require(MASS)
#
## START
if(verbose)cat("\n## HAPLIN, VERSION 3.0.1 ##\n")
#
## READ DATA AND REPORT MISSING:
if(verbose)	cat("\nReading data from file...  ")
	if((design == "triad") | (design == "cc.triad")) {
		.fam <- "mfc"
	}
	if(design == "cc") .fam <- "c"
	.data.read <- f.read.data(indata = .info$filename, sep = .info$filespecs$sep, allele.sep = .info$filespecs$allele.sep, na.strings = .info$filespecs$na.strings, markers = .info$filespecs$markers, use.missing = .info$model$use.missing, variables = .info$filespecs$n.vars, family = .fam) ##
	
###	return(.data.read)
	
if(verbose)	cat("Done\n")

	.rows.with.na <- attr(.data.read, "rows.with.na")
	.rows.dropped <- attr(.data.read, "rows.dropped")

#
	.ntri.seq <- rep(NA, 4)
	names(.ntri.seq) <- c("Original", "After rem NA", "After rem Mend. inc.", "After rem unused haplos")
	.ntri.seq[2] <- dim(.data.read)[1]
	
	if(.rows.with.na == 0){
		if(verbose) cat("No lines contained missing data\n")	
		.ntri.seq[1] <- .ntri.seq[2]
	}
	else {
		if(use.missing){
			if(verbose) cat("There were ", .rows.with.na, " rows with missing data\nAll rows retained in analysis\n", sep = "")
			.ntri.seq[1] <- .ntri.seq[2]
			}
		else{
			if(verbose) cat("The following", .rows.with.na, "data lines were dropped due to missing data:\n", .rows.dropped, "\n")
			.ntri.seq[1] <- .ntri.seq[2] + .rows.with.na
		}
	}

#
## FREQUENCY COUNT AND ALLELE SORTING:
if(verbose) cat("\nPreparing data for analysis...  ")
##	.data0 <- f.prep.data(.data.read)

	.data <- f.prep.data.new(.data.read, short = T)	

	
if(verbose) cat("Done\n")

## EXTRACT ALLELE INFORMATION:
.info$haplos$alleles <- attr(.data, "alleles")
#
## CHANGE CASE: UPPER-CASE IS MOST FREQUENT	
.f.change.case <- function(allele){
	names(allele) <- casefold(names(allele), upper = F)
	.max <- which(allele == max(allele))[1]
	names(allele)[.max] <- casefold(names(allele)[.max], upper = T)
	allele
}
.info$haplos$alleles <- lapply(.info$haplos$alleles, .f.change.case) # IS THIS A GOOD IDEA?
#	
##
.nloci <- length(.info$haplos$alleles)
#
#


## DESIGN-DEPENDENT DATA PREPARATIONS:
##
#

.tmp <- f.sep.data(.data, .info)
.data.gen <- .tmp$data.gen
.data.vars <- .tmp$data.vars
.HWE.res <- .tmp$HWE.res




.orig.lines <- attr(.data.gen, "orig.lines") # A LIST OF THE ORIGINAL LINE NUMBERS. CAN BE INDEXED BY ind.unique.line


## CHECK SOME OF THE HWE RESULTS:
for(i in seq(along = .info$haplos$alleles)) if(any(.info$haplos$alleles[[i]]  != .HWE.res[[i]]$freq)) warning("Something's strange with the frequency count in HWE test!")


if(design == "triad" | design == "cc.triad"){
#
## REPORT MENDELIAN INCONSISTENCIES:	
	.rows.with.Mendelian.inconsistency <- attr(.data.gen, "rows.with.Mendelian.inconsistency")
#
	if(length(.rows.with.Mendelian.inconsistency) == 0){
		.ind.Mend <- numeric(0)
		.ntri.seq[3] <- .ntri.seq[2]
		if(!use.missing & .rows.with.na > 0)
			if(verbose) cat("None of the retained lines contained Mendelian inconsistencies\n")
		else
			if(verbose) cat("No lines contained Mendelian inconsistencies\n")
	}
	else
	{
		if(use.missing | .rows.with.na == 0) .ind.Mend <- .rows.with.Mendelian.inconsistency
		else{
		.ind.Mend <- seq(length = dim(.data.read)[1] + .rows.with.na)
		.ind.Mend <- .ind.Mend[-.rows.dropped][.rows.with.Mendelian.inconsistency]
		}
		if(verbose) cat("The following", length(.ind.Mend), "data lines were dropped due to Mendelian inconsistencies:\n", .ind.Mend, "\n")	
		.ntri.seq[3] <- .ntri.seq[2] - length(.ind.Mend)
	}
}
if(design == "cc"){
	.ntri.seq[3] <- .ntri.seq[2] # CC CANNOT DETECT MEND. INCONS...
}
#
#
##
## PRELIMINARY DATA FIXUP:
#
#	EXPAND FREQUENCIES AND ADD COUNTER:
		.orig.lines <- unlist(.orig.lines[.data.gen$ind.unique.line])
		
		.ind <- 1:(dim(.data.gen)[1])
		.ind <- rep(.ind, .data.gen$freq)
		.ind.aux <- unlist(sapply(.data.gen$freq, function(x)1:x))
#
##
	if(design == "triad" | design == "cc.triad"){
		if(!xchrom){
			.data.gen <- cbind(.data.gen[.ind,1:5], ind.aux = .ind.aux, .orig.lines)
			names(.data.gen) <- c("m1", "m2", "f1", "f2", "ind.unique.line", "ind.aux", "orig.lines")
		}
		if(xchrom){
			.data.gen <- cbind(.data.gen[.ind,1:6], ind.aux = .ind.aux, .orig.lines)
			names(.data.gen) <- c("m1", "m2", "f1", "f2", "sex", "ind.unique.line", "ind.aux", "orig.lines")
		}

##		row.names(.data.gen) <- NULL
	}
##
	if(design == "cc"){
###		.data.gen <- cbind(.data.gen[.ind,1:5], ind.aux = .ind.aux)
		.data.gen <- cbind(.data.gen[.ind,1:3], ind.aux = .ind.aux, .orig.lines)
###		names(.data.gen) <- c("m1", "m2", "f1", "f2", "ind.unique.line", "ind.aux")
		names(.data.gen) <- c("c1", "c2", "ind.unique.line", "ind.aux", "orig.lines")
##		row.names(.data.gen) <- NULL
	}

##	REPLACE LINE COUNTERS ETC. WITH UNIQUE TAG ind, WHICH HAS ONE VALUE FOR 
##	EACH (REMAINING) TRIAD:
.tag.tmp <- f.create.tag(.data.gen[,c("ind.unique.line", "ind.aux")])
.tag.tmp <- match(.tag.tmp, unique(.tag.tmp))
.data.gen$ind <- .tag.tmp
.data.gen$ind.unique.line <- .data.gen$ind.aux <- NULL
###attr(.data.gen, "alleles") <- .info$haplos$alleles # DENNE ANTAGELIG IKKE LENGRE NODVENDIG
#



.data <- .data.gen # TEMPORARY!!

if(xchrom) return(.data)


#
## COMPUTE PRELIMINARY HAPLOTYPE FREQUENCIES USING A SIMPLE EM-VERSION:
##
.prelim.freq <- f.preliminary.freq.new(.data, .info)
.data$freq <- .prelim.freq
.info$haplos$prelim.haplotype.freq <- attr(.prelim.freq, "prelim.haplotype.freq")
#
## DECIDE WHICH HAPLOTYPES TO INCLUDE IN ANALYSIS
.info$haplos$selected.haplotypes <- f.sel.haplos(.info)
.n.sel.haplos <- sum(.info$haplos$selected.haplotypes)

#
#
## REMOVE OR JOIN (NO, THAT OPTION HAS BEEN REMOVED...) HAPLOTYPES WITH INITIAL FREQUENCY BELOW threshold.
## HAPLOTYPES ARE RECODED TO 1,2,3,... AFTER REMOVAL.
## FREQUENCIES ARE RENORMALIZED SO THAT EACH TRIAD SUMS TO ONE.
##
if(verbose) cat("Removing unused haplotypes...  ")
	
if(abs(sum(.data$freq) - .ntri.seq[3]) > 1e-6) warning("There may be a problem with the data summary")
.data <- f.repl.thin(.data, selection = .info$haplos$selected.haplotypes, design = design)
attr(.data, "selected.haplotypes") <- .info$haplos$selected.haplotypes # BURDE IKKE VÆRE NØDVENDIG....
.ntri.seq[4] <- sum(.data$freq)
if(abs(.ntri.seq[4] - round(.ntri.seq[4])) > 1e-6) warning("There may be a problem with the data summary")
.ntri.seq[4] <- round(.ntri.seq[4])
if(verbose) cat("Done\n")

#
## DECIDE REFERENCE
.tmp <- f.prep.reference(.info)
reference.method <- .tmp$reference.method
ref.cat <- .tmp$ref.cat
.info$haplos$ref.cat <- ref.cat

#
#
## ADD ON CASE-CONTROL VARIABLE FOR cc AND cc.triad DATA
if(design == "cc" | design == "cc.triad"){
	.ccvar <- .info$variables$ccvar
	.cc <- .data.vars[.data$orig.lines, .ccvar]
	if(any(is.na(.cc))) stop(paste(sum(is.na(.cc)), " missing values found in case-control variable! Must be removed from file before analysis!\n", sep = ""))
	.codes <- names(attr(.data.vars, "variables")[[.ccvar]])
	if(length(.codes) != 2) stop(paste('Case-control variable "ccvar" is coded with ', paste(.codes, collapse = ", "), '. It should have exactly two different values!', sep = ""))	
	if(!identical(sort(unique(.cc)), c(1,2))) stop("Something's wrong with the case-control variable!") # SHOULDN'T BE NECESS.
	if(verbose) cat("\nNote: The following case/control coding has been assumed:\ncontrols = ", .codes[1], ", cases = ", .codes[2], "\n", sep = "")
	# if(!identical(.codes, c("0","1"))) stop("Case-control variable must be coded as 0 (control) and 1 (case)!")
	# .cc <- as.numeric(.codes[.cc])
	.data$cc <- .cc
}
#
#
## START ESTIMATION, INCLUDING EM IF REQUIRED AND RESAMPLING IF REQUESTED:
.use.EM <- T # PROBABLY NO REASON NOT TO USE EM
if(!.use.EM){
	stop("not yet checked!\n")
	.agg.data <- f.sum.and.expand(.data)
	.res.0 <- f.tri.glm(.agg.data$freq, maternal = maternal, ref.cat = ref.cat, response = "simple", design = design)		
	.res <- f.tri.glm(.agg.data$freq, maternal = maternal, ref.cat = ref.cat, response = .info$haplos$response, design = design)
}
if(.use.EM){
	if(verbose) cat("\nUsing EM to estimate model with no effect:\n")
	.res.0 <- f.EM.missing(data = .data, maternal = maternal, verbose = verbose, ref.cat = ref.cat, design = design, response = "simple", max.EM.iter = max.EM.iter)
	.info$estimation$iter.used$iter.used.0 <- attr(.res.0, "iter.used")
	.info$estimation$EM.conv$EM.conv.0 <- attr(.res.0, "EM.conv")
	if(verbose)cat("\nDone\n")
#
	if(verbose) cat("\nUsing EM to estimate full model:\n")
#
	if(.info$model$scoretest == "only" & !.info$control$data.out){
		## RUN ONLY ONE ITERATION, TO OBTAIN X
		## MERK: EGENTLIG NOK BARE Å FÅ UT X. 
		## DESSUTEN: BURDE HER KUNNE KJØRE BÅDE mult OG free. MEN: free GÅR IKKE NØDVENDIGVIS FOR cc
		.res.oneiter <- f.EM.missing(data = .data, response = .info$haplos$response, maternal = maternal, verbose = verbose, ref.cat = ref.cat, design = design, max.EM.iter = 0, x = T, suppressEmWarnings = T)
	}else{
		.res <- f.EM.missing(data = .data, response = .info$haplos$response, maternal = maternal, verbose = verbose, ref.cat = ref.cat, design = design, max.EM.iter = max.EM.iter, x = T)
		.info$estimation$iter.used$iter.used <- attr(.res, "iter.used")
		.info$estimation$EM.conv$EM.conv <- attr(.res, "EM.conv")
	}
	if(verbose)cat("\nDone\n")
} # END if(.use.EM)
#
.info$estimation$iter.used <- c(iter.used.0 = attr(.res.0, "iter.used"), iter.used = attr(.res, "iter.used"))
.info$estimation$EM.conv <- c(EM.conv.0 = attr(.res.0, "EM.conv"), EM.conv = attr(.res, "EM.conv"))
#
if(resampling == "jackknife"){
	if(design != "triad") stop("Jackknifing has not been tested with case-control data....")
	if(.info$model$scoretest == "only") stop('Jackknifing has not been tested when scoretest == "only"')
	if(verbose) cat("\nStarting jackknife\n")
	.res.resamp <- f.jackknife.new(data = .data, maternal = maternal, ref.cat = ref.cat, verbose = F, use.EM = .use.EM, max.EM.iter = max.EM.iter)
	attr(.res, "cov.resamp") <- .res.resamp$cov
}
#
## DATA OUT, IF REQUESTED
if(data.out){
	if(design == "triad" | design == "cc.triad"){
		.pos <- f.pos.in.grid(A = rep(.n.sel.haplos, 4), comb = as.matrix(.data[,c("m1", "m2", "f1", "f2")]))
		if(design == "cc.triad") cat("Er dette virkelig rett for cc.triad? (ser litt rart ut i farten)\n")
	}
	if(design == "cc"){
		.pos <- f.pos.in.grid(A = c(rep(.n.sel.haplos, 2), 2), comb = as.matrix(.data[,c("c1", "c2", "cc")]))
	}
		.pred <- .res$pred[.pos]
		.predsum <- f.groupsum(X = .pred, INDICES = .data$ind)
		.freqsum <- f.groupsum(X = .data$freq, INDICES = .data$ind) # MERK: .freqsum ER 1 FOR DENNE VARIANTEN
		.pred.redist <- .pred/.predsum * .freqsum 
	# (KAN HENDE DENNE ER TILGJENGELIG SOM y I OBJECT!!)

	return(cbind(.data, pred = .pred.redist))
}


if(F & design != "triad") {
## 
	f.vis(.HWE.res, vis = T)
	f.vis(.prelim.haplotype.freq, vis = T)
	.tmput <- .prelim.haplotype.freq * 2 * sum(.data$freq)
	names(.tmput) <- .haplotypes
	
	f.vis(.tmput, vis = T)

	f.vis(.info$haplos$alleles, vis = T)
print(summary(.res, reference.method = reference.method, design = "cc"), design = "cc")
	
	return(.res)
###	return(.prelim.haplotype.freq * 2 * sum(.data$freq))
	return(.data)
}
#
if(.info$model$scoretest %in% c("yes", "only")){
	 ## COMPUTE VAR-COVAR AND SCORE UNDER NULL HYPOTHESIS:
	if(.info$model$scoretest == "yes") .restemp <- .res
	if(.info$model$scoretest == "only") .restemp <- .res.oneiter
	.tmp.0 <- f.var.covar(pred = .res.0$pred, X = .restemp$result$x, data = .data, info = .info)
	#
	.score.0 <- .tmp.0[["score"]]
	.var.cov.0 <- .tmp.0[["var.covar"]]
	.nullpars <- length(.res.0$result$coefficients) # KAN VEL FINNE BEDRE MÅTE?
	.score.0.red <- .score.0[,-(1:.nullpars), drop = F]
	.var.cov.0.red <- .var.cov.0[-(1:.nullpars), -(1:.nullpars), drop = F]
	#
	## SUM SCORE
	.sc.0 <- t(.score.0.red) %*% rep(1, dim(.score.0.red)[1]) # SUM SCORE OVER INDIVIDUAL FAMILIES
	#
	## CHI-SQUARED VALUE & TEST
	.sc.test <- as.numeric(t(.sc.0) %*% .var.cov.0.red %*% .sc.0)
	.sc.df <- length(.restemp$result$coefficients) - .nullpars ## HAR VEL DF'S ALLEREDE?
	.sc.pval <- pchisq(.sc.test, df = .sc.df, lower.tail = F) 
	#cat("##############\n")
	#print(.sc.pval)
	#cat("##############\n")
	#
	.score.ut <- list(score = .score.0, chisquared = .sc.test, df = .sc.df, pval = .sc.pval)
	#
	if(.info$model$scoretest == "only"){
		## DUMP INFORMATION
		.score.ut <- c(list(info = .info), .score.ut)
		return(.score.ut)
	}
}else{
	.score.ut <- NULL
	.var.cov.0 <- NULL
}
#
## JUST CHECK THAT THE TRIAD ACCOUNTING IS CORRECT
if(max(abs(.res$ntri - .ntri.seq[4]), abs(.res.0$ntri - .ntri.seq[4])) > 1e-6) warning("There may be a problem with the data summary")
#
## COMPUTE "EXACT" VAR-COVAR, TAKING EM INTO ACCOUNT
.var.cov <- f.var.covar(pred = .res$pred, X = .res$result$x, data = .data, info = .info)[["var.covar"]]
attr(.res, "cov.correct") <- .var.cov

#
#f.vis(round(diag(.var.cov), 4))
#f.vis(round(diag(summary(.res$result)$cov.unscaled), 4))
#
## OBTAIN LOG-LIKELIHOOD, VIA DIFFERENT METHODS (SHOULD BE SORTED OUT...):
.anova <- anova(.res.0$result, .res$result, test = "Chisq")
.df <- .anova$Df[2]

if(F){
	.lratio.test <- anova(.res.0$result, .res$result, test = "Chisq")[2,"P(>|Chi|)"]
}	
f.vis(.loglike.ratio.full <- .anova$Deviance[2], vis = F)

.loglike.0 <- f.final.loglike(data = .data, pred = .res.0$pred, info = .info)
.loglike <- f.final.loglike(data = .data, pred = .res$pred, info = .info)

f.vis(.loglike.ratio <- 2*(.loglike - .loglike.0), vis = F)

.loglike.ratio <- c(.loglike.ratio, full.loglike = .loglike.ratio.full)
names(.loglike.ratio) <- paste(names(.loglike.ratio), ".ratio", sep = "")

.lratio.test <- 1 - pchisq(.loglike.ratio, df = .df)

.lratio.test <- .lratio.test[1] ## BRUKER BARE DEN FORSTE, DE ANDRE VAR EKSPERIMENTELLE

.lratio.ut <- c(loglike.0 = unname(.loglike.0[1]), loglike = unname(.loglike[1]), df = .df, p.value.overall = unname(.lratio.test))
	
## return(.lratio.test)	

if(F){	
	f.vis(.tmpan <- anova(.res.0$result, .res$result, test = "Chisq"), vis = T)
	print(anova(.res.0$result, .res$result, test = "Chisq"), digits = 12)

	# f.vis(.temp.0, vis = T)
	# f.vis(.temp, vis = T)

	# f.vis(2*(.temp[1] - .temp.0[1]), vis = T)
	# print(2*(.temp[1] - .temp.0[1]), digits = 12)

	# f.vis(2*(.temp[2] - .temp.0[2]), vis = T)
	# print(2*(.temp[2] - .temp.0[2]), digits = 12)

	## stop()

	f.vis(2*(.temp - .temp.0), vis = T)
	print(2*(.temp - .temp.0), digits = 12)

	f.vis(.res.0$result$dev - .res$result$dev, vis = T)
}
#
if(verbose) cat("\nEstimation finished, preparing output...  ")
.var.cov.ut <- list(var.cov.0 = .var.cov.0, var.cov = .var.cov)
#
#
#.out <- list(result = .res, design = design, alleles = .info$haplos$alleles, selected.haplotypes = .info$haplos$selected.haplotypes, resampling = resampling, orig.call = sys.call(), date = date(), reference.method = reference.method, rows.dropped = .rows.dropped, HWE.res = .HWE.res, ntri.seq = .ntri.seq, loglike = .lratio.ut, score = .score, temp = .var.cov.0, temp.multi = .tmp.ut$var.covar.multinom, info = .info)
#.out <- list(result = .res, design = design, alleles = .info$haplos$alleles, selected.haplotypes = .info$haplos$selected.haplotypes, resampling = resampling, orig.call = sys.call(), date = date(), reference.method = reference.method, rows.dropped = .rows.dropped, HWE.res = .HWE.res, ntri.seq = .ntri.seq, loglike = .lratio.ut, score = .score.ut, temp = .var.cov.0, info = .info)
.out <- list(result = .res, design = design, alleles = .info$haplos$alleles, selected.haplotypes = .info$haplos$selected.haplotypes, resampling = resampling, orig.call = sys.call(), date = date(), reference.method = reference.method, rows.dropped = .rows.dropped, HWE.res = .HWE.res, ntri.seq = .ntri.seq, loglike = .lratio.ut, score = .score.ut, var.cov = .var.cov.ut, info = .info)
class(.out) <- "haplin"	#
if(verbose) cat("Done\n")
#
if(printout){
## PRINTOUT AND PLOTTING, FOR THE CONVENIENCE OF THE USER:
	plot(.out, reference = reference.method)#
#
	if(verbose) cat("\n#################################\n")
	#
	print(summary(.out, reference = reference.method)) #
}
	invisible(.out)
}
