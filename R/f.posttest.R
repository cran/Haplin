"f.posttest"<-
function(object.list, test = c("single", "double"))
{
# TEST FOR DIFFERENCE IN PARAMETER ESTIMATES OVER A SERIES OF HAPLIN OBJECTS
# WARNING: THE SET OF HAPLOTYPES USED IN EACH ESTIMATION *MUST* BE THE SAME!
# object.list IS A LIST OF HAPLIN OBJECTS. 
# test CAN INCLUDE ANY OF "haplo.freq", "single", "double"
#### PREPARE: ###################
require(MASS)
#
## TEST FOR CONSISTENCY BETWEEN ELEMENTS OF LIST
cat("BOER SJEKKE MER FOR KONSISTENS, F.EKS. AT maternal = F FOR ALLE!\n")
if(length(object.list) <= 1) stop("Need at least two elements in object.list")
.info <- lapply(object.list, function(x) x$info)
#
.selected.haplotypes <- lapply(.info, function(x)x$haplos$selected.haplotypes)
.sel.haps <- .selected.haplotypes[[1]]
.sjekk <- sapply(.selected.haplotypes[-1], function(x) identical(x, .sel.haps))
if(any(!.sjekk)) stop()
#
.ref.cats <- lapply(.info, function(x)x$haplos$ref.cat)
.ref.cat <- .ref.cats[[1]]
.sjekk <- sapply(.ref.cats[-1], function(x) identical(x, .ref.cat))
if(any(!.sjekk)) stop()
#
## EXTRACT SEPARATE RESULTS AND COEFFICIENTS
.res <- lapply(object.list, function(x) x$result$result)
.coef <- lapply(.res, function(x) x$coefficients)
.names <- names(.coef[[1]])
#
.f.cov <- function(res){
## EXTRACT COVARIANCE MATRIX
## PASS PAA! HVILKEN MATRISE BLIR BRUKT?!
	.summary <- summary.glm(res, correlation = F)
	.cov <- .summary$cov.unscaled
	#
	.cov.correct <- attr(res, "cov.correct")
	if(!is.null(.cov.correct)) .cov <- .cov.correct
	#
	.cov.resamp <- attr(res, "cov.resamp")
	if(!is.null(.cov.resamp)) .cov <- .cov.resamp
	return(.cov)
}
#
.cov <- lapply(.res, .f.cov)
#
## FIND POSITIONS OF RELEVANT PARAMETERS
.mf <- grep("mf", .names, value = T)
.c <- grep("c[[:digit:]]", .names, value = T)
.cdd <- grep("cdd[[:digit:]]", .names, value = T)
if((length(.mf) == 0) | (length(.c) == 0) | (length(.cdd) == 0))stop("Something's wrong with the coefficient names")
#
## FOR THE HAPLOTYPE FREQUENCIES, SUBTRACT FIRST PARAMETER FROM THE REST,
## TO "NORMALIZE".
.tmp <- f.post.diff(coeff = .coef, covar = .cov)
.coef <- .tmp$coeff
.cov <- .tmp$cov
#
.names <- .names[-1]
.mf <- .mf[-1]
#
## EXTRACT RELEVANT PARAMETERS
.coef.c <- lapply(.coef, function(x)x[.c, , drop = F])
.coef.cdd <- lapply(.coef, function(x)x[.cdd, , drop = F])
.coef.comb <- lapply(.coef, function(x)x[c(.c, .cdd), , drop = F])
#
.cov.c <- lapply(.cov, function(x) x[.c, .c, drop = F])
.cov.cdd <- lapply(.cov, function(x) x[.cdd, .cdd, drop = F])
.cov.comb <- lapply(.cov, function(x) x[c(.c, .cdd), c(.c, .cdd), drop = F])
#
.contr.c <- diag(length(.coef.c[[1]]))
.contr.cdd <- diag(length(.coef.cdd[[1]]))
.contr.comb <- diag(length(.coef.comb[[1]]))
#
.res.c <- vector(length(.coef.c), mode = "list")
.res.cdd <- vector(length(.coef.cdd), mode = "list")
.res.comb <- vector(length(.coef.comb), mode = "list")

for (i in seq(along = .coef.c)){
	.res.c[[i]] <- f.post.chisq(coeff = .coef.c[[i]], covar = .cov.c[[i]], contrast.mat = .contr.c)
	.res.cdd[[i]] <- f.post.chisq(coeff = .coef.cdd[[i]], covar = .cov.cdd[[i]], contrast.mat = .contr.cdd)
	.res.comb[[i]] <- f.post.chisq(coeff = .coef.comb[[i]], covar = .cov.comb[[i]], contrast.mat = .contr.comb)
	f.vis(.res.c, vis = F)
	f.vis(.res.cdd, vis = F)
}
cat("\nIndividual Wald tests within each stratum, \nfor single dose, double dose and combined:\n")
cat("\nPost-hoc (Wald) test single dose:")
.res.c.vis <- cbind("Stratum: ", seq(along = .res.c), ", Chi-squared = ", round(sapply(.res.c, function(x)x$chisq), 3), ", df's = ", sapply(.res.c, function(x) x$df), ", p-value = ", round(sapply(.res.c, function(x) x$pval), 5))
dimnames(.res.c.vis) <- list(rep("", dim(.res.c.vis)[1]), rep("", dim(.res.c.vis)[2]))
print(.res.c.vis, quote = F, print.gap = 0)

cat("\nPost-hoc (Wald) test double dose:")
.res.cdd.vis <- cbind("Stratum: ", seq(along = .res.cdd), ", Chi-squared = ", round(sapply(.res.cdd, function(x)x$chisq), 3), ", df's = ", sapply(.res.cdd, function(x) x$df), ", p-value = ", round(sapply(.res.cdd, function(x) x$pval), 5))
dimnames(.res.cdd.vis) <- list(rep("", dim(.res.cdd.vis)[1]), rep("", dim(.res.cdd.vis)[2]))
print(.res.cdd.vis, quote = F, print.gap = 0)

cat("\nPost-hoc (Wald) test combined single and double dose:")
.res.comb.vis <- cbind("Stratum: ", seq(along = .res.comb), ", Chi-squared = ", round(sapply(.res.comb, function(x)x$chisq), 3), ", df's = ", sapply(.res.comb, function(x) x$df), ", p-value = ", round(sapply(.res.comb, function(x) x$pval), 5))
dimnames(.res.comb.vis) <- list(rep("", dim(.res.comb.vis)[1]), rep("", dim(.res.comb.vis)[2]))
print(.res.comb.vis, quote = F, print.gap = 0)

cat("\nCompare combined to overall likelihood ratio p-values:")
.p.value.overall <- sapply(object.list, function(x){
	.tmp <- summary(x)$loglike["p.value.overall"]
	if(is.null(.tmp)) .tmp <- NA
	return(.tmp)
})
.p.value.overall.vis <- cbind("Stratum: ", seq(along = .p.value.overall), ", p-value = ", round(as.numeric(.p.value.overall), 5))
dimnames(.p.value.overall.vis) <- list(rep("", nrow(.p.value.overall.vis)), rep("", ncol(.p.value.overall.vis)))
print(.p.value.overall.vis, quote = F, print.gap = 0)


#####################
#
## SELECT PARAMETERS TO BE TESTED
.names.list <- list(haplo.freq = .mf, single = .c, double = .cdd)
.nam <- names(.names.list)
if(!all(test %in% .nam)) stop('Invalid input in argument "test"')
f.vis(.velg <- unlist(.names.list[.nam %in% test])) # MAKE SURE SELECTION IS IN CORRECT ORDER
#
f.vis(.coef <- lapply(.coef, function(x)x[.velg, , drop = F]))
f.vis(.cov <- lapply(.cov, function(x) x[.velg, .velg, drop = F]))
#
## RESHAPE COEFFICIENTS AND COVARIANCE MATR. INTO FULL SIZE
.n.pars <- length(.coef[[1]])
.l <- length(.coef)
f.vis(.coef.vec <- unlist(.coef))
f.vis(.cov.mat <- f.bdiag(.cov))
#
## BUILD CONTRAST MATRIX
.A <- f.post.contrasts(test.type = "interaction", n.res = .l, n.pars = .n.pars)
#
## DO CHI-SQUARED TEST
.chisq.res <- f.post.chisq(coeff = .coef.vec, covar = .cov.mat, contrast.mat = .A)

cat("\nWald test of heterogeneity\n")
cat("Tested effects: '", paste(test, collapse = "' '"), "'", sep = "")
.Wald.vis <- cbind(c("Chi-squared value:", "Df's:", "P-value:"), round(c(.chisq.res$chisq, .chisq.res$df, .chisq.res$pval), 5))
dimnames(.Wald.vis) <- list(rep("", dim(.Wald.vis)[1]), rep("", dim(.Wald.vis)[2]))
print(.Wald.vis, quote = F, print.gap = 0)

#print(round(.chisq.res, 5))



return(invisible(.chisq.res))

return(.chisq)
return(.A)
return(.cov)
return()
f.vis(.coef, .cov, vis = T)







.n.haplos <- sum(.sel.haps)

### BURDE SUNNHETSSJEKKE .valg!


print(.coef)
print(.cov)

### HVORDAN VAR DET MED INVERTERING, SA DU??

.inv <- lapply(.cov, function(x) solve(x))
.chisq <- rep(NA, length(.coef))

for(i in seq(along = .coef)){
	.chisq[i] <- .coef[[i]] %*% .inv[[i]] %*% .coef[[i]]
}
.df <- length(.coef) * length(.coef[[1]])
.chisq <- sum(.chisq)
.pval <- pchisq(.chisq, df = .df, lower.tail = F)

print(.pval)


return(.chisq)
return(.inv)



return()



.coef.alt <- lapply(.res, function(x) summary(x)$coefficients)

print(.coef)
print(.coef.alt)


return(.coef)

}
