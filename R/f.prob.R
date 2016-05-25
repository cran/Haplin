f.prob <- function(grid.sim, haplo.freq, RRcm, RRcf, RRstar, sim.maternal, RR.mat, RRstar.mat, sim.xchrom, sim.comb.sex = "double", BR.girls){
##
## grid.sim: GRID OF HAPLOTYPE (AND SEX) COMBINATIONS
## haplo.freq: HAPLOTYPE FREQUENCIES
## RRcm, RRcf: VECTOR OF RELATIVE RISKS, MATERNALLY AND PATERNALLY INHERITED
## RRstar IS A VECTOR OF R* IN MODEL
## sim.maternal, ARGUMENT TO INDICATE IF MATERNAL EFFECTS SHOULD BE INHERITED
## RR.mat AND RRstar.mat MATERNAL RR AND RRstar
## sim.xchrom: INDICATES WHETHER SIMULATION IS ON xchrom
## sim.comb.sex: SIMILAR TO comb.sex IN haplin
## BR.girls: BASELINE RISK FOR GIRLS (RELATIVE TO BOYS)
##
#
if(sim.xchrom) .sex <- grid.sim[, "sex"]
#
## HAPLOTYPE FREQUENCIES
.hf1 <- haplo.freq[grid.sim[,1]]
.hf2 <- haplo.freq[grid.sim[,2]]
.hf3 <- haplo.freq[grid.sim[,3]]
if(!sim.xchrom) .hf4 <- haplo.freq[grid.sim[,4]]
#
## RELATIVE RISKS, SINGLE DOSE
.RRcm <- RRcm[grid.sim[,"h2.m"]]
.RRcf <- RRcf[grid.sim[,"h2.f"]]	
#
## MATERNAL CONTRIBUTIONS
if(sim.maternal){
	.RR1.mat <- RR.mat[grid.sim[,"h1.m"]]
	.RR2.mat <- RR.mat[grid.sim[,"h2.m"]]
}
#
## RELATIVE RISKS, DOUBLE DOSE
.RRstar <- rep(1, length(.RRcm)) # default 1
 # set to RRstar where same allele from mother and father
.ind.double <- grid.sim[,"h2.m"] == grid.sim[,"h2.f"]
.RRstar[.ind.double] <- RRstar[grid.sim[,"h2.m"]][.ind.double]
#
## XCHROM ADAPTATIONS
if(sim.xchrom){
	## GIRLS ARE UNCHANGED, BOYS GET MODIFIED
	if(sim.comb.sex == "single"| sim.comb.sex == "females"| sim.comb.sex == "males"){
		## NO CONTRIBUTION FROM FATHERS TO SONS
		.RRcf[.sex == 1] <- 1 
		.RRstar[.sex == 1] <- 1
	}
	if(sim.comb.sex == "double"){
		## FOR BOYS, MAKE SURE CONTRIBUTION "FROM FATHER" MATCHES THAT FROM MOTHER
		## SO THAT EFFECT OF A SINGLE ALLELE FROM MOTHER IS EQUAL TO RRcm * RRcf * RRstar
		.RRcf[.sex == 1] <- RRcf[grid.sim[,"h2.m"]][.sex == 1]
		.RRstar[.sex == 1] <- RRstar[grid.sim[,"h2.m"]][.sex == 1]
	}
}
if(sim.maternal){
	.RRstar.mat <- rep(1, length(.RR1.mat))
	.ind.double.mat <- grid.sim[,"h1.m"] == grid.sim[,"h2.m"]
	.RRstar.mat[.ind.double.mat] <- RRstar.mat[grid.sim[,"h2.m"]][.ind.double.mat]
}
#
## MULTINOMIAL PROBABILITIES, ASSUMING HWE AND MULTIPLICATIVE RISKS
if(sim.xchrom){
	## SEPARATE BASELINE LEVEL FOR GIRLS, AND ONLY THREE ALLELES, IF xchrom
	if(sim.comb.sex=="females") .prob <- .hf1 * .hf2 * .hf3 * .RRcm * .RRcf * .RRstar * (.sex == 2)
	else if (sim.comb.sex=="males") .prob <- .hf1 * .hf2 * .hf3 * .RRcm * .RRcf * .RRstar * (.sex == 1)
	else .prob <- .hf1 * .hf2 * .hf3 * .RRcm * .RRcf * .RRstar * ((.sex == 1) + (.sex == 2) * BR.girls)
}else{
	.prob <- .hf1 * .hf2 * .hf3 * .hf4 * .RRcm * .RRcf * .RRstar
}
if(sim.maternal){
	.prob <- .prob * .RR1.mat * .RR2.mat * .RRstar.mat
}
#
## NORMALIZE
.prob <- .prob/sum(.prob)
#
return(.prob)
#
}
