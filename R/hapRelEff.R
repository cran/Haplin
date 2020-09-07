hapRelEff <- function(nall = 2, cases.comp, controls.comp, cases.ref, controls.ref, haplo.freq, RR, RRcm, RRcf, RRstar, RR.mat, RRstar.mat, xchrom = F, sim.comb.sex = "double", BR.girls, response = "mult", ...){	
#
## Compute relative efficiency for two different study designs 
#
## Mir: Indep of alpha, that fits with Pitman-eff., but not with Bahadur... Not relevant to compare different designs with different levels, always?
#
## Number of possible haplotypes
.nhaplo <- prod(nall)
#
## The relative efficiency is only calculated for one stratum
n.strata <- 1
##
##
## Make cases and controls list arguments
## Unnecessary to keep the arguments as lists when resticted to one stratum, but makes it easier to expand to several strata later on.
## Also in agreement with hapPowerAsymp etc...
if(missing(controls.comp)) controls.comp <- list(c(mfc=0))
if(missing(controls.ref)) controls.ref <- list(c(mfc=0))
if(!is.list(cases.comp)) cases.comp <- list(cases.comp)
if(!is.list(cases.ref)) cases.ref <- list(cases.ref)
if(!is.list(controls.comp)) controls.comp <- list(controls.comp)
if(!is.list(controls.ref)) controls.ref <- list(controls.ref)
#
## Misc tests
#
if(response=="free") stop("response = \"free\" is not yet implemented", call. = F)
#
.sim.maternal <- FALSE
if(!missing(RR.mat) | !missing(RRstar.mat)) .sim.maternal <- TRUE
.sim.poo <- FALSE
if(!missing(RRcm) | !missing(RRcf)) .sim.poo <- TRUE
if(.sim.poo && !missing(RR)) stop("RR cannot be present at the same time as RRcm and RRcf", call. = F)
#
## Not in use
if(F){
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
}
#
## 
if(missing(RRstar)) RRstar <- rep(1,.nhaplo)
if(.sim.maternal && missing(RRstar.mat)) RRstar.mat <- rep(1,.nhaplo)
#
## Arguments to be passed on to hapCovar
.asymp.arg <- c(as.list(environment()), list(...))
.asymp.arg1 <- .asymp.arg2 <- .asymp.arg
.asymp.arg1[which(names(.asymp.arg1)%in%c("cases.ref","controls.ref"))] <- NULL
.asymp.arg2[which(names(.asymp.arg2)%in%c("cases.comp","controls.comp"))] <- NULL
names(.asymp.arg1)[which(names(.asymp.arg1)=="cases.comp")] <- "cases"
names(.asymp.arg1)[which(names(.asymp.arg1)=="controls.comp")] <- "controls"
names(.asymp.arg2)[which(names(.asymp.arg2)=="cases.ref")] <- "cases"
names(.asymp.arg2)[which(names(.asymp.arg2)=="controls.ref")] <- "controls"
#
## Find asymptotic covariance matrices from hapCovar
.asymp1 <- do.call(hapCovar, args=.asymp.arg1)
.asymp2 <- do.call(hapCovar, args=.asymp.arg2)
#
.cov1 <- .asymp1$cov
.cov2 <- .asymp2$cov
.coef <- .asymp1$coef #.asymp1$coef is equal to .asymp2$coef
#
.ref.cat <- attr(.asymp1,"ref.cat")
#
## NAMES OF ALL COEFFICIENTS
.names <- rownames(.coef[[1]])
## SPLIT IN EFFECT GROUPS
.effs <- f.coefnames(.names)
#
## NAMES OF COEFFICIENTS ASSOCIATED WITH EACH TEST
.names.tests <- list(haplo.freq = .effs$haplo.freq, child = c(.effs$child.s, .effs$child.d), child.poo = c(.effs$child.poo.m, .effs$child.poo.f, .effs$child.d),  poo = .effs$poo, maternal = c(.effs$maternal.s, .effs$maternal.d))
#
.covar1 <- lapply(unique(unlist(.names.tests)), function(x){
	f.vis(.coef <- lapply(.coef, function(y) y[x, , drop = F]), vis = F)
	f.vis(.coef.vec <- unlist(.coef), vis = F)
	f.vis(.cov <- lapply(.cov1, function(y) y[x, x, drop = F]), vis = F)
	f.vis(.cov.mat <- f.bdiag(.cov), vis = F)
	.contrast.mat = diag(1,length(.coef[[1]]))
	.cov <- .contrast.mat %*% .cov.mat %*% t(.contrast.mat)
})
.covar2 <- lapply(unique(unlist(.names.tests)), function(x){
	f.vis(.coef <- lapply(.coef, function(y) y[x, , drop = F]), vis = F)
	f.vis(.coef.vec <- unlist(.coef), vis = F)
	f.vis(.cov <- lapply(.cov2, function(y) y[x, x, drop = F]), vis = F)
	f.vis(.cov.mat <- f.bdiag(.cov), vis = F)
	.contrast.mat = diag(1,length(.coef[[1]]))
	.cov <- .contrast.mat %*% .cov.mat %*% t(.contrast.mat)
})
#
## Number of genotyped individuals
.n1 <- cases.comp[[1]]*nchar(names(cases.comp[[1]])) + controls.comp[[1]]*nchar(names(controls.comp[[1]]))
.n2 <- cases.ref[[1]]*nchar(names(cases.ref[[1]])) + controls.ref[[1]]*nchar(names(controls.ref[[1]]))
#.n1 <- sum(sapply(1:n.strata, function(x) cases.comp[[x]]*nchar(names(cases.comp[[x]])) + controls.comp[[x]]*nchar(names(controls.comp[[x]]))))
#.n2 <- sum(sapply(1:n.strata, function(x) cases.ref[[x]]*nchar(names(cases.ref[[x]])) + controls.ref[[x]]*nchar(names(controls.ref[[x]]))))
#
## Relative efficieny
.rel.eff <- sapply(1:length(unique(unlist(.names.tests))), function(x){
	rel.eff <- (.covar2[[x]]*.n2)/(.covar1[[x]]*.n1)
	return(rel.eff)
})
#
names(.rel.eff) <- unique(unlist(.names.tests))
#
## Mir: Leftovers from overall-test. Not sure if output should be a list
.rel.eff <- list(.rel.eff)
names(.rel.eff) <- "haplo.rel.eff"

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
	.out.haplo$RRcm_cf.rel.eff <- .out.haplo$RRcf.rel.eff <- .out.haplo$RRcm.rel.eff <- rep("ref",.nhaplo)
	.out.haplo$RRcm.rel.eff[-.ref.cat] <- round(.rel.eff[["haplo.rel.eff"]][which(grepl("cm",names(.rel.eff[["haplo.rel.eff"]])) & !grepl("_",names(.rel.eff[["haplo.rel.eff"]]))), drop=F],2)
	.out.haplo$RRcf.rel.eff[-.ref.cat] <- round(.rel.eff[["haplo.rel.eff"]][which(grepl("cf",names(.rel.eff[["haplo.rel.eff"]])) & !grepl("_",names(.rel.eff[["haplo.rel.eff"]]))), drop=F],2)
	.out.haplo$RRcm_cf.rel.eff[-.ref.cat] <- round(.rel.eff[["haplo.rel.eff"]][which(grepl("_",names(.rel.eff[["haplo.rel.eff"]]))), drop=F],2)
} else{
	.out.haplo$RR.rel.eff <- rep("ref",.nhaplo)
	.out.haplo$RR.rel.eff[-.ref.cat] <- round(.rel.eff[["haplo.rel.eff"]][which(grepl("c",names(.rel.eff[["haplo.rel.eff"]])) & !grepl("dd",names(.rel.eff[["haplo.rel.eff"]]))), drop=F],2)
}
if(.sim.maternal){
	.out.haplo$RRm.rel.eff <- rep("ref",.nhaplo)
	.out.haplo$RRm.rel.eff[-.ref.cat] <- round(.rel.eff[["haplo.rel.eff"]][which(!grepl("c",names(.rel.eff[["haplo.rel.eff"]])) & !grepl("dd",names(.rel.eff[["haplo.rel.eff"]])) & !grepl("mf",names(.rel.eff[["haplo.rel.eff"]]))), drop=F],2)	
}
#
.out <- list(.out.haplo)
names(.out) <- "haplo.rel.eff"
#
return(.out)
#
}

