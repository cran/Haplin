hapPowerAsymp <- function(nall=2, n.strata = 1, cases, controls, haplo.freq, RR, RRcm, RRcf, RRstar, RR.mat, RRstar.mat, xchrom = F, sim.comb.sex = "double", BR.girls, response = "mult", alpha = 0.05, ...){
#
#
if(response=="free") stop("response = \"free\" is not yet implemented", call. = F)
#
## Number of possible haplotypes
.nhaplo <- prod(nall)
#
if(alpha<0 | alpha>1) stop("The significance level, alpha, must have a value between 0 and 1", call. = F)
#
.sim.maternal <- FALSE
if(!missing(RR.mat) | !missing(RRstar.mat)) .sim.maternal <- TRUE
.sim.poo <- FALSE
if(!missing(RRcm) | !missing(RRcf)) .sim.poo <- TRUE
#
if(missing(RRstar)) RRstar <- rep(1,.nhaplo)
if(.sim.maternal && missing(RRstar.mat)) RRstar.mat <- rep(1,.nhaplo)
#
## Save misc tests for hapCovar
#
## Arguments to be passed on to hapCovar
.asymp.arg <- c(as.list(environment()), list(...))
#
## Find asymptotic covariance matrix from hapCovar
.asymp <- do.call(hapCovar, args=.asymp.arg)
.cov <- .asymp$cov
.coef <- .asymp$coef
.ref.cat <- attr(.asymp,"ref.cat")
#
## NAMES OF ALL COEFFICIENTS
.names <- rownames(.coef[[1]])
## SPLIT IN EFFECT GROUPS
.effs <- f.coefnames(.names)
#
## Choose tests to be run
.test <- "haplo.freq"
if(.sim.poo){ # POO
	.child.poo <- F # works, but probably not relevant for standard use
	#
	if(.child.poo) .test <- c(.test, "child.poo")
	.test <- c(.test, "poo")
} else { # NOT POO
	.test <- c(.test, "child")
}
if(.sim.maternal) .test <- c(.test, "maternal")
#
## NAMES OF COEFFICIENTS ASSOCIATED WITH EACH TEST
#if(response=="free" & !.sim.poo) .names.tests <- list(haplo.freq = .effs$haplo.freq, child = c(.effs$child.s, .effs$child.d), maternal = c(.effs$maternal.s, .effs$maternal.d))
#else if(response=="free" & .sim.poo) .names.tests <- list(haplo.freq = .effs$haplo.freq, child.poo = c(.effs$child.poo.m, .effs$child.poo.f, .effs$child.d),  poo = .effs$poo, maternal = c(.effs$maternal.s, .effs$maternal.d))
.names.tests <- list(haplo.freq = .effs$haplo.freq, child = c(.effs$child.s, .effs$child.d), child.poo = c(.effs$child.poo.m, .effs$child.poo.f, .effs$child.d),  poo = .effs$poo, maternal = c(.effs$maternal.s, .effs$maternal.d))
#
## Find lambda, i.e., the chi-squared value
if(n.strata!=1){
	.chisq.res <- lapply(unique(unlist(.names.tests)), function(x) f.posttest(coef_ = .coef, cov_ = .cov, test = x))
} else	.chisq.res <- lapply(unique(unlist(.names.tests)), function(x){
	## SELECT COEFFICIENTS TO BE TESTED
	f.vis(.coef <- lapply(.coef, function(y) y[x, , drop = F]), vis = F)
	f.vis(.cov <- lapply(.cov, function(y) y[x, x, drop = F]), vis = F)
	## RESHAPE COEFFICIENTS AND COVARIANCE MATR. INTO FULL SIZE
	f.vis(.coef.vec <- unlist(.coef), vis = F)
	f.vis(.cov.mat <- f.bdiag(.cov), vis = F)
	## DO CHI-SQUARED TEST
	.chisq.res <- f.post.chisq(coeff = .coef.vec, covar = .cov.mat, contrast.mat = diag(1,length(.coef[[1]])))
})
#
names(.chisq.res) <- unique(unlist(.names.tests))
#
## Remove redundant y-vector and reduce to numeric
.chisq.res <- lapply(.chisq.res, function(x){
	x$y <- NULL
	return(unlist(x))
})
#
## If .nhaplo > 2, find overall p-value
if(.nhaplo>2){
	if(n.strata!=1){
		.chisq.res.overall <- lapply(.test, function(x) f.posttest(coef_ = .coef, cov_ = .cov, test = .names.tests[[x]]))
	} else	.chisq.res.overall <- lapply(.test, function(x){
		## SELECT COEFFICIENTS TO BE TESTED
		f.vis(.coef <- lapply(.coef, function(y) y[.names.tests[[x]], , drop = F]), vis = F)
		f.vis(.cov <- lapply(.cov, function(y) y[.names.tests[[x]], .names.tests[[x]], drop = F]), vis = F)
		## RESHAPE COEFFICIENTS AND COVARIANCE MATR. INTO FULL SIZE
		f.vis(.coef.vec <- unlist(.coef), vis = F)
		f.vis(.cov.mat <- f.bdiag(.cov), vis = F)
		## DO CHI-SQUARED TEST
		.chisq.res.overall <- f.post.chisq(coeff = .coef.vec, covar = .cov.mat, contrast.mat = diag(1,length(.coef[[1]])))
	})
	#
	names(.chisq.res.overall) <- .test
	#
	## Remove redundant y-vector and reduce to numeric
	.chisq.res.overall <- lapply(.chisq.res.overall, function(x){
		x$y <- NULL
		return(unlist(x))
	})
}
#
## Calculate power from non-centrality parameter 
.power <- sapply(.chisq.res, function(x){1-pchisq(qchisq((1-alpha),x["df"], ncp=0), x["df"], ncp = x["chisq"])})
if(.nhaplo>2) .power.overall <- sapply(.chisq.res.overall, function(x){1-pchisq(qchisq((1-alpha),x[["df"]], ncp=0), x[["df"]], ncp = x[["chisq"]])})
#
## Prepare output
#
## Haplotypes
.haplo <- lapply(nall,function(x) 1:x)
.haplo <- do.call("expand.grid",.haplo)
.names <- names(.haplo)
if(length(nall)>1) .haplo <- apply(.haplo[,.names],1,paste,collapse="-")
.out.haplo <- data.frame(.haplo, stringsAsFactors=F)
names(.out.haplo) <- "Haplotype"
#
if(.sim.poo){
	.out.haplo$RRcm_cf.power <- .out.haplo$RRcf.power <- .out.haplo$RRcm.power <- rep("ref",.nhaplo)
	.out.haplo$RRcm.power[-.ref.cat] <- round(.power[which(grepl("cm",names(.power)) & !grepl("_",names(.power))), drop=F],2)
	.out.haplo$RRcf.power[-.ref.cat] <- round(.power[which(grepl("cf",names(.power)) & !grepl("_",names(.power))), drop=F],2)
	.out.haplo$RRcm_cf.power[-.ref.cat] <- round(.power[which(grepl("_",names(.power))), drop=F],2)
} else{
	.out.haplo$RR.power <- rep("ref",.nhaplo)
	.out.haplo$RR.power[-.ref.cat] <- round(.power[which(grepl("c",names(.power)) & !grepl("dd",names(.power))), drop=F],2)
}
if(.sim.maternal){
	.out.haplo$RRm.power <- rep("ref",.nhaplo)
	.out.haplo$RRm.power[-.ref.cat] <- round(.power[which(!grepl("c",names(.power)) & !grepl("dd",names(.power)) & !grepl("mf",names(.power))), drop=F],2)	
}
#
.out <- list(.out.haplo)
names(.out) <- "haplo.power"
if(.nhaplo>2){
	.power.overall <- .power.overall[-1]
	.out <- list(.out.haplo,.power.overall)
	names(.out) <- c("haplo.power","overall.power")
}
#
return(.out)
}