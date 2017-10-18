f.hapSim <- function(nall, cases, controls, haplo.freq, RR, RRcm, RRcf, RRstar, sim.maternal, RR.mat, RRstar.mat, sim.poo, gen.missing.cases, gen.missing.controls, xchrom, sim.comb.sex, BR.girls, nhaplo, nloci, case.des, control.des, case.mat, control.mat, case.design, control.design){
##
## To be used by function hapSim and hapRun
## Simulates genetic data in Haplin format, consisting of fetal effects, maternal effects and/or parent-of-origin effects. 
## Allows for simulation of both autosomal and X-linked markers, assuming HWE and multiplicative risks.
##
##
#
#
if(sim.poo){
	.RRcm <- RRcm
	.RRcf <- RRcf
}else{
	.RRcm <- .RRcf <- RR
}
#
## Generate missing values at random
gen.missing <- function(gen.missing, .file, .l){
	if(length(gen.missing) == 1) .gen.missing <- rep(gen.missing,((ncol(.file)-.l+1)/2))
	else{ 
		if(is.vector(gen.missing)) .gen.missing <- rep(gen.missing,each = 3)
		else .gen.missing <- as.vector(t(gen.missing))
	}
	.j <- 1 
	for(i in seq(.l,ncol(.file),by=2)){
		.row <- which(rbinom(nrow(.file),1,.gen.missing[.j]) == 1)
		.j <- .j+1
		.file[.row,c(i,i+1)] <- NA
	}
	return(.file)
}	
#
## Relative risks for controls
.RR.controls <- rep(1,nhaplo)
#
## Basic grid of haplotypes:
if(xchrom){
	.grid <- as.matrix(expand.grid(h1.m = 1:nhaplo, h2.m = 1:nhaplo, h2.f = 1:nhaplo, sex = c(1,2)))
}else{
	.grid <- as.matrix(expand.grid(h1.m = 1:nhaplo, h2.m = 1:nhaplo, h1.f = 1:nhaplo, h2.f = 1:nhaplo))
}
#
## Probability calculations
.prob.cases <- f.prob(grid.sim = .grid, haplo.freq = haplo.freq, RRcm = .RRcm, RRcf = .RRcf, RRstar = RRstar, sim.maternal = sim.maternal, RR.mat = RR.mat, RRstar.mat = RRstar.mat, sim.xchrom = xchrom, sim.comb.sex = sim.comb.sex, BR.girls = BR.girls)
.prob.controls <- f.prob(grid.sim = .grid, haplo.freq = haplo.freq, RRcm = .RR.controls, RRcf = .RR.controls, RRstar = .RR.controls, sim.maternal = FALSE, RR.mat = .RR.controls, RRstar.mat = .RR.controls, sim.xchrom = xchrom, sim.comb.sex = sim.comb.sex, BR.girls = BR.girls)
#	
## Simulate cases
.file.cases <- as.matrix(f.sim(.prob = .prob.cases, size = rowSums(case.mat), nall = nall, .nloci = nloci, xchrom = xchrom, .grid = .grid))
#
## Set columns missing (if family members are missing by design)
## The number of columns to the left of the genetic data equals .l-1
.l <- 1
if(xchrom) .l <- 2
.k <- 1
.n <- case.mat[1]
for(j in 1:ncol(case.mat)){ 
	if(!grepl("m", case.design[j]))	.file.cases[.k:.n, c(seq(.l, 6*nloci, 6),seq(.l+1, 6*nloci, 6))] <- NA
	if(!grepl("f", case.design[j]))	.file.cases[.k:.n, c(seq(.l+2, 6*nloci, 6),seq(.l+3, 6*nloci, 6))] <- NA
	.k <- .n + 1
	.n <- .n + case.mat[j+1]
}
#
## Genereate missing values at random
if(!is.null(gen.missing.cases)) .file.cases <- gen.missing(gen.missing = gen.missing.cases, .file = .file.cases, .l = .l)
.file.cases <- cbind(rep(1, rowSums(case.mat)), .file.cases)
colnames(.file.cases)[1] <- "cc"
#		
## Simulate controls and bind controls to cases 
if(!rowSums(control.mat) == 0){
	#
	## Set columns missing (if family members are missing by design)
	.k <- 1
	.n <- control.mat[1]
	.file.controls <- as.matrix(f.sim(.prob = .prob.controls, size = rowSums(control.mat), nall = nall, .nloci = nloci, xchrom = xchrom, .grid = .grid))
	for(j in 1:ncol(control.mat)){
		if(!grepl("m", control.design[j])) .file.controls[.k:.n, c(seq(.l, 6*nloci, 6), seq(.l+1, 6*nloci, 6))] <- NA	
		if(!grepl("f", control.design[j])) .file.controls[.k:.n, c(seq(.l+2, 6*nloci, 6), seq(.l+3, 6*nloci, 6))] <- NA
		if(!grepl("c", control.design[j])) .file.controls[.k:.n, c(seq(.l+4, 6*nloci, 6), seq(.l+5, 6*nloci, 6))] <- NA	
		.k <- .n + 1
		.n <- .n + control.mat[j+1]
	}
	#
	## Genereate missing values at random
	if(!(is.null(gen.missing.controls) | rowSums(control.mat) == 0)) .file.controls <- gen.missing(gen.missing = gen.missing.controls, .file = .file.controls, .l)
	#
	.file.controls <- cbind(rep(0, rowSums(control.mat)), .file.controls)
	colnames(.file.controls)[1] <- "cc"
	#
	## Bind cases and controls by rows
	.file <- rbind(.file.cases, .file.controls)
	#
} else .file <- .file.cases
#
## Return invisible
return(invisible(as.data.frame(.file, stringsAsFactors=FALSE)))
}
