f.hapArg <- function(nall, n.strata, cases, controls, haplo.freq, RR, RRcm, RRcf, RRstar, RR.mat, RRstar.mat, gen.missing.cases=NULL, gen.missing.controls=NULL, n.sim, xchrom=F, sim.comb.sex, BR.girls, sim.maternal, sim.poo, nloci, nhaplo){
##
##
#
.n.strata <- n.strata
#
.arg <- as.list(match.call())
.arg[1] <- NULL
.arg$n.strata <- NULL
#
if(is.null(gen.missing.cases) | is.matrix(gen.missing.cases)) .arg$gen.missing.cases <- list(gen.missing.cases)
if(is.null(gen.missing.controls) | is.matrix(gen.missing.controls)) .arg$gen.missing.controls <- list(gen.missing.controls)
#
.strat.arg <- lapply(.arg, function(x){
	if(!is.list(x) & is.vector(x)) .strat.arg <- rep(list(x),.n.strata)
	else if(length(x) == 1) .strat.arg <- rep(x,.n.strata)
	else if(length(x) == .n.strata) .strat.arg <- x
	else .strat.arg <- NA	
	}
)
#
## Misc tests
if(is.list(cases) & !is.null(names(cases))) stop("Argument cases is not specified correctly", call. = F)
if(is.list(controls) & !is.null(names(controls))) stop("Argument controls is not specified correctly", call. = F)
if(any(sapply(.strat.arg,function(x){length(x)==1 && is.na(x)}))) stop(paste("The argument",names(which(sapply(.strat.arg,function(x){length(x)==1 && is.na(x)}))),"is not specified correctly \n", sep = " "), call. = F)
#
## Possible family designs
.case.des <- c("mfc", "mc", "fc", "c")
.control.des <- c("mfc", "mc", "fc", "mf", "c", "m", "f")
#
.cc.arg <- lapply(1:.n.strata, function(x){
	.case.mat <- t(as.matrix(.strat.arg$cases[[x]]))
	.control.mat <- t(as.matrix(.strat.arg$controls[[x]]))
	.case.design <- names(.strat.arg$cases[[x]])
	.control.design <- names(.strat.arg$controls[[x]])
	.cc.arg <- list(case.des=.case.des, control.des=.control.des, case.mat=.case.mat, control.mat=.control.mat, case.design=.case.design, control.design=.control.design)
})
#
.strat.arg <- cbind(do.call("cbind",.strat.arg), do.call("rbind",.cc.arg))
#
return(.strat.arg)
#
}
