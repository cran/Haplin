f.hapTests <- function(nall, cases, controls, case.design, control.design, case.mat, control.mat, case.des, control.des, haplo.freq, RR, RRcm, RRcf, RRstar, sim.maternal, sim.poo, RR.mat, RRstar.mat, gen.missing.cases, gen.missing.controls, xchrom, sim.comb.sex, BR.girls, nloci, nhaplo, n.sim){
##
## Misc tests to be run by hapSim and hapRun
##
#
if(n.sim < 1 | n.sim%%1 != 0) stop("Argument \"nsim\" must be an integer equal to or larger than 1", call. = F)
#
if(any(nall < 2) | any(nall%%1 != 0)) stop("Argument \"nall\" must be integer(s) larger than 1", call. = F)
#
if(!(sim.comb.sex %in% c("double","single","females","males"))) stop("\"sim.comb.sex\" is not specified correctly", call. = F)
#
if(!all(case.design %in% case.des) || any(duplicated(names(cases)))) stop("Argument \"cases\" is not specified correctly", call. = F)
if(!is.numeric(case.mat) || any(case.mat <= 0) || any(case.mat%%1 != 0)) stop("Argument \"cases\" must contain positive integers", call. = F)
if(!all(control.design %in% control.des) || any(duplicated(names(controls)))) stop("Argument \"controls\" is not specified correctly", call. = F)
if(!is.numeric(control.mat) || any(control.mat < 0) || any(control.mat%%1 != 0)) stop("Argument \"controls\" must contain only non-negative integers", call. = F)
#
if(length(control.mat) > 1 && any(control.mat==0)) stop("Argument \"controls\" is not specified correctly", call. = F)
#
if(sim.poo){
	if(!all(c(length(haplo.freq),length(RRcm),length(RRcf),length(RRstar)) == nhaplo)) stop("The relative risks and haplotype frequencies must have length equal to the number of haplotypes", call. = F)
	if(any(c(RRcm,RRcf,RRstar)<0)) stop("The relative risks cannot have negative values", call. = F)
}else{
	if(!all(c(length(haplo.freq),length(RR),length(RRstar)) == nhaplo)) stop("The relative risks and haplotype frequencies must have length equal to the number of haplotypes", call. = F)
	if(any(c(RR,RRstar)<0)) stop("The relative risks cannot have negative values", call. = F)
}	
if(sim.maternal){
	if(!all(c(length(haplo.freq), length(RR.mat), length(RRstar.mat)) == nhaplo)) stop("The relative risks and haplotype frequencies must have length equal to the number of haplotypes", call. = F)
	if(any(c(RR.mat,RRstar.mat)< 0)) stop("The relative risks cannot have negative values", call. = F)
}
#
if(xchrom && BR.girls <= 0) stop("Argument \"BR.girls\" must be positive", call. = F)
#
if(!is.null(gen.missing.cases) && !(is.numeric(gen.missing.cases) && all(gen.missing.cases <= 1 & gen.missing.cases >= 0 & (length(gen.missing.cases) == 1 || length(gen.missing.cases) == nloci || (is.matrix(gen.missing.cases) && all(dim(gen.missing.cases) == c(nloci, 3))))))) stop("Argument \"gen.missing.cases\" is not specified correctly", call. = F)
if(!is.null(gen.missing.controls) & rowSums(control.mat) == 0) warning("There are no control families present and the argument \"gen.missing.controls\" is thus ignored", call. = F)
else if(!is.null(gen.missing.controls) && !(is.numeric(gen.missing.controls) && all(gen.missing.controls <= 1 & gen.missing.controls >= 0 & (length(gen.missing.controls) == 1 || length(gen.missing.controls) == nloci || (is.matrix(gen.missing.controls) && all(dim(gen.missing.controls) == c(nloci, 3))))))) stop("Argument \"gen.missing.controls\" is not specified correctly", call. = F)
if(!all(grepl("m",case.design) & grepl("f",case.design) & grepl("c",case.design)) && length(gen.missing.cases) == 3*nloci) warning("The values in \"gen.missing.cases\" corresponding to family members that are missing by design are ignored", call. = F)
if(!all(grepl("m",control.design) & grepl("f",control.design) & grepl("c",control.design)) && length(gen.missing.controls) == 3*nloci) warning("The values in \"gen.missing.controls\" corresponding to family members that are missing by design are ignored", call. = F)
#
return(invisible())
#
}
