postTest <- function(object.list)
{
# TEST FOR DIFFERENCE IN PARAMETER ESTIMATES OVER A SERIES OF HAPLIN OBJECTS
# WARNING: THE SET OF HAPLOTYPES USED IN EACH ESTIMATION _MUST_ BE THE SAME!
# object.list IS A LIST OF HAPLIN OBJECTS. 
# test.haplo CAN INCLUDE ANY OF "haplo.freq", "single", "double"
#
# test = c("child.single"), 
#test = c("c", "cdd")
#### PREPARE: ###################
.vis <- F
#
## PRELIMINARY CHECKS
if(any(is.na(object.list))) return(object.list)# JUST PASS RIGHT THROUGH IF AT LEAST ONE NA 
if(class(object.list) != "haplinStrat") stop("Argument 'object.list' should be the result from running 'haplinStrat'", call. = F)
#
## REMOVE OVERALL RESULT, ONLY DO TESTING ON SUBSTRATA, OF COURSE
.object.list <- object.list[-1]
if(length(.object.list) <= 1) stop("Need at least two strata", call. = F)
#
## STRATA NAMES, TO BE USED IN PRINTOUT
.stratnavn <- names(.object.list)
if(is.null(.stratnavn)) .stratnavn <- as.character(seq(along = .object.list)) ## TRENGS IKKE
#
## EXTRACT MISC INFO
.info <- lapply(.object.list, function(x) x$info)
.response <- .info[[1]]$haplos$response
.maternal <- .info[[1]]$model$maternal
.poo <- .info[[1]]$model$poo

#
## CONSISTENCY CHECK OF HAPLOTYPES AMONG ELEMENTS OF LIST
.tmp.selected.haplotypes <- lapply(.info, function(x)x$haplos$selected.haplotypes)
.sel.haps <- .tmp.selected.haplotypes[[1]]
.sjekk <- sapply(.tmp.selected.haplotypes[-1], function(x) identical(tolower(x), tolower(.sel.haps)))
if(any(!.sjekk)) stop("Different haplotypes selected in different strata!",
call. = F)
#
## CONSISTENCY CHECK OF ref.cat AMONG ELEMENTS OF LIST
.tmp.ref.cat <- lapply(.info, function(x)x$haplos$ref.cat)
.ref.cat <- .tmp.ref.cat[[1]]
.sjekk <- sapply(.tmp.ref.cat[-1], function(x) identical(x, .ref.cat))
if(any(!.sjekk)) stop()

#
## EXTRACT SEPARATE RESULTS, COEFFICIENTS, AND COVAR-MATRICES
.params <- lapply(.object.list, coef)
.coef <- lapply(.params, function(x) x$coef)
.cov <- lapply(.params, function(x) x$cov)
#
## FOR THE HAPLOTYPE FREQUENCIES, SUBTRACT FIRST PARAMETER FROM THE REST,
## TO "NORMALIZE". 
.tmp <- f.post.diff(coeff = .coef, covar = .cov)
.coef <- .tmp$coeff
.cov <- .tmp$cov
#
.names <- rownames(.coef[[1]])
#
## FIND NAMES/POSITIONS OF RELEVANT PARAMETERS. NOTE: \\< INSISTS ON START OF WORD, SO THAT, FOR INSTANCE, "cm1" ISN'T PICKED UP BY "m"
.mf <- grep("\\<mf", .names, value = T)
.c <- grep("\\<c[[:digit:]]", .names, value = T)
.cm <- grep("\\<cm[[:digit:]]", .names, value = T)
.cf <- grep("\\<cf[[:digit:]]", .names, value = T)
.cdd <- grep("\\<cdd[[:digit:]]", .names, value = T)
.m <- grep("\\<m[[:digit:]]", .names, value = T)
.mdd <- grep("\\<mdd[[:digit:]]", .names, value = T)
# SOME AD HOC TESTING
.flag <- (length(.mf) == 0) |
	(!.poo & (length(.c) == 0)) |
	(.poo & ((length(.cm) == 0) | (length(.cf) == 0))) |
	((length(.cdd) == 0) & (.response == "free")) |
	(.maternal && (length(.m) == 0)) |
	(.maternal && (.response == "free") && (length(.mdd) == 0))
if(.flag) stop("Something's wrong with the coefficient names", call. = F)
#
## 
standard.tests <- F
if(standard.tests){
	#
	## EXTRACT RELEVANT PARAMETERS
	.f.extr <- function(co, selpars){
		if(ncol(co[[1]]) == 1){
			## COEFFICIENTS
			.res <- lapply(co, function(x) x[selpars, , drop = F])
		}else {
			## COVARIANCE MAT
			.res <- lapply(co, function(x) x[selpars, selpars, drop = F])
		}
			return(.res)
	}
	.coef.c <- .f.extr(.coef, .c)
	if(.response == "free") .coef.cdd <- .f.extr(.coef, .cdd)
	.coef.comb <- .f.extr(.coef, c(.c, .cdd))
	#
	.cov.c <- .f.extr(.cov, .c)
	if(.response == "free") .cov.cdd <- .f.extr(.cov, .cdd)
	.cov.comb <- .f.extr(.cov, c(.c, .cdd))
	#
	## CONTRAST MATRICES
	.contr.c <- diag(length(.coef.c[[1]]))
	if(.response == "free") .contr.cdd <- diag(length(.coef.cdd[[1]]))
	.contr.comb <- diag(length(.coef.comb[[1]]))
	#
	## RESULT LISTS
	.res.c <- vector(length(.coef.c), mode = "list")
	if(.response == "free") .res.cdd <- vector(length(.coef.cdd), mode = "list")
	.res.comb <- vector(length(.coef.comb), mode = "list")

	for (i in seq(along = .coef.c)){
		.res.c[[i]] <- f.post.chisq(coeff = .coef.c[[i]], covar = .cov.c[[i]], contrast.mat = .contr.c)
		if(.response == "free") .res.cdd[[i]] <- f.post.chisq(coeff = .coef.cdd[[i]], covar = .cov.cdd[[i]], contrast.mat = .contr.cdd)
		.res.comb[[i]] <- f.post.chisq(coeff = .coef.comb[[i]], covar = .cov.comb[[i]], contrast.mat = .contr.comb)
		f.vis(.res.c, vis = F)
		if(.response == "free") f.vis(.res.cdd, vis = F)
	}
	cat("\nIndividual Wald tests within each stratum, \nfor single dose, double dose and combined:\n")
	cat("\nPost-hoc (Wald) test single dose:")
	.res.c.vis <- cbind("Stratum: ", format(.stratnavn, justify = "right"), ", Chi-squared = ", round(sapply(.res.c, function(x)x$chisq), 3), ", df's = ", sapply(.res.c, function(x) x$df), ", p-value = ", round(sapply(.res.c, function(x) x$pval), 5))
	dimnames(.res.c.vis) <- list(rep("", dim(.res.c.vis)[1]), rep("", dim(.res.c.vis)[2]))
	print(.res.c.vis, quote = F, print.gap = 0)

	if(.response == "free") {
		cat("\nPost-hoc (Wald) test double dose:")
		.res.cdd.vis <- cbind("Stratum: ", format(.stratnavn, justify = "right"), ", Chi-squared = ", round(sapply(.res.cdd, function(x)x$chisq), 3), ", df's = ", sapply(.res.cdd, function(x) x$df), ", p-value = ", round(sapply(.res.cdd, function(x) x$pval), 5))
		dimnames(.res.cdd.vis) <- list(rep("", dim(.res.cdd.vis)[1]), rep("", dim(.res.cdd.vis)[2]))
		print(.res.cdd.vis, quote = F, print.gap = 0)
	}
	cat("\nPost-hoc (Wald) test combined single and double dose:")
	.res.comb.vis <- cbind("Stratum: ", format(.stratnavn, justify = "right"), ", Chi-squared = ", round(sapply(.res.comb, function(x)x$chisq), 3), ", df's = ", sapply(.res.comb, function(x) x$df), ", p-value = ", round(sapply(.res.comb, function(x) x$pval), 5))
	dimnames(.res.comb.vis) <- list(rep("", dim(.res.comb.vis)[1]), rep("", dim(.res.comb.vis)[2]))
	print(.res.comb.vis, quote = F, print.gap = 0)

	cat("\nCompare combined to overall likelihood ratio p-values:")
	.p.value.overall <- sapply(.object.list, function(x){
		.tmp <- summary(x)$loglike["p.value.overall"]
		if(is.null(.tmp)) .tmp <- NA
		return(.tmp)
	})
	.p.value.overall.vis <- cbind("Stratum: ", format(.stratnavn, justify = "right"), ", p-value = ", round(as.numeric(.p.value.overall), 5))
	dimnames(.p.value.overall.vis) <- list(rep("", nrow(.p.value.overall.vis)), rep("", ncol(.p.value.overall.vis)))
	print(.p.value.overall.vis, quote = F, print.gap = 0)
}

#####################
#
## DO THE INTERACTION TESTING FOR HAPLO FREQUENCIES, CHILD EFFECTS AND, IF RELEVANT, MATERNAL AND CHILD+MATERNAL
.chisq.res <- vector(2, mode = "list")
names(.chisq.res) <- c("haplo.freq", "child")
#
.chisq.res[["haplo.freq"]] <- f.posttest(coef_ = .coef, cov_ = .cov, mf = .mf, c_ = .c, cm_ = .cm, cf = .cf, cdd = .cdd, m = .m, mdd = .mdd, test = "haplo.freq")
if(.poo){
	.chisq.res[["poo"]] <- f.posttest(coef_ = .coef, cov_ = .cov, mf = .mf, c_ = .c, cm_ = .cm, cf = .cf, cdd = .cdd, m = .m, mdd = .mdd, test = "poo")
	#
	if(.maternal){
		.chisq.res[["maternal"]] <- f.posttest(coef_ = .coef, cov_ = .cov, mf = .mf, c_ = .c, cm_ = .cm, cf = .cf, cdd = .cdd, m = .m, mdd = .mdd, test = "maternal")
		.chisq.res[["poo.and.mat"]] <- f.posttest(coef_ = .coef, cov_ = .cov, mf = .mf, c_ = .c, cm_ = .cm, cf = .cf, cdd = .cdd, m = .m, mdd = .mdd, test = c("poo", "maternal"))
	}
}else{
	.chisq.res[["child"]] <- f.posttest(coef_ = .coef, cov_ = .cov, mf = .mf, c_ = .c, cm_ = .cm, cf = .cf, cdd = .cdd, m = .m, mdd = .mdd, test = "child")
	#
	if(.maternal){
		.chisq.res[["maternal"]] <- f.posttest(coef_ = .coef, cov_ = .cov, mf = .mf, c_ = .c, cm_ = .cm, cf = .cf, cdd = .cdd, m = .m, mdd = .mdd, test = "maternal")
		.chisq.res[["chi.and.mat"]] <- f.posttest(coef_ = .coef, cov_ = .cov, mf = .mf, c_ = .c, cm_ = .cm, cf = .cf, cdd = .cdd, m = .m, mdd = .mdd, test = c("child", "maternal"))
	}
}



.ut <- lapply(.chisq.res, function(x) {
	x$y <- NULL
	return(unlist(x))
})

.ut <- do.call("rbind", .ut)
.ut <- dframe(test = rownames(.ut), .ut)
rownames(.ut) <- NULL

return(.ut)
#
#	return(.chisq.res)
#
#
#	cat("\nWald test of heterogeneity\n")
#	cat("Tested effects: '", paste(test, collapse = "' '"), "'", sep = "")
#	.Wald.vis <- cbind(c("Chi-squared value:", "Df's:", "P-value:"), round(c(.chisq.res$chisq, .chisq.res$df, .chisq.res$pval), 5))
#	dimnames(.Wald.vis) <- list(rep("", dim(.Wald.vis)[1]), rep("", dim(.Wald.vis)[2]))
#	print(.Wald.vis, quote = F, print.gap = 0)


}

