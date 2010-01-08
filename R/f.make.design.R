f.make.design <- function(

#observed, maternal = F, ref.cat, response = "free", design = "triad", xchrom = F, ...
n.alleles, maternal, ref.cat, response, design, xchrom,

test = F

)
{
#
# THE PROGRAM ESTIMATES EFFECTS OF SEVERAL ALLELES IN A CASE-TRIAD, CASE-CONTROL-TRIAD OR CASE-CONTROL DESIGN
# ASSUMING HARDY-WEINBERG EQUILIBRIUM (AND RARE DISEASE, IF NECESSARY)
# IMPORTANT: OBSERVED FREQUENCIES (observed) MUST COME FROM A COMPLETE
# GRID IN APPROPRIATE ORDER! THE ORDER IS:
# TRIAD: BY MOTHER'S FIRST ALLELE AND SECOND ALLELE,
# THEN THE FATHER'S FIRST AND SECOND
# CASE-CONTROL-TRIAD: BY MOTHER'S FIRST ALLELE AND SECOND ALLELE,
# THEN THE FATHER'S FIRST AND SECOND, AND FINALLY CONTROL (1) AND CASE (2), 
# SO THAT TOP HALF OF THE DATA ARE CONTROLS
# CASE-CONTROL: BY FIRST ALLELE, SECOND ALLELE
# AND FINALLY CONTROL (1) AND CASE (2), SO THAT TOP HALF OF THE DATA ARE
# CONTROLS
#
# maternal DETERMINES WHETHER OR NOT MATERNAL EFFECTS SHOULD BE ESTIMATED
# (NOT AVAILABLE FOR THE CASE-CONTROL DESIGN)
# THE ... PASSES ON OTHER ARGUMENTS, LIKE start, TO THE GLM


#
# CREATE GRID OF PARENTAL (OR FETAL) GENOTYPES (AND/OR POSSIBLY CASE-CONTROL STATUS)
.A <- n.alleles
.tmp <- factor(1:.A)
if((design == "triad") & !xchrom){
	.mf <- expand.grid(m1 = .tmp, m2 = .tmp, f1 = .tmp, f2 = .tmp)	#
}
if((design == "triad") & xchrom){
	.mf <- expand.grid(m1 = .tmp, m2 = .tmp, f2 = .tmp, sex = c(0,1)) ## HERE, 0 AND 1 (MALES AND FEMALES, RESP.) ARE USED AS DUMMY. ORIGINAL CODING (OUTSIDE f.tri.glm) SHOULD BE 1 (MALES) AND 2 (FEMALES)
}
if(design == "cc.triad"){	
	if(xchrom)stop("Not implemented")
	.mf <- expand.grid(m1 = .tmp, m2 = .tmp, f1 = .tmp, f2 = .tmp, cc = c(0,1)) ## HERE, 0 AND 1 (CONTROL AND CASE, RESP.) ARE USED AS DUMMY. ORIGINAL CODING (OUTSIDE f.tri.glm) SHOULD BE 1 (CONTROL) AND 2 (CASE). IN INPUT FILE, THE LARGEST IS ALWAYS THE CASE
}
if(design == "cc"){	
	if(xchrom)stop("Not implemented")
	.mf <- expand.grid(m2 = .tmp, f2 = .tmp, cc = c(0,1)) ## HERE, 0 AND 1 (CONTROL AND CASE, RESP.) ARE USED AS DUMMY. ORIGINAL CODING (OUTSIDE f.tri.glm) SHOULD BE 1 (CONTROL) AND 2 (CASE)
}
#
## CREATE DUMMY VARIABLES FOR PATERNAL GENOTYPES
if(design != "cc"){
	.m1.dum <- model.matrix( ~ -1 + m1, data = .mf)
	if(!xchrom) .f1.dum <- model.matrix( ~ -1 + f1, data = .mf)
}
.m2.dum <- model.matrix( ~ -1 + m2, data = .mf) # DUMMIES FOR MATERNALLY DERIVED HAPLOTYPE
.f2.dum <- model.matrix( ~ -1 + f2, data = .mf) # DUMMIES FOR PATERNALLY DERIVED HAPLOTYPE
#
## CREATE "DUMMY" (COUNTING) VARIABLES FOR CHILDREN:
.mm.c <- .m2.dum + .f2.dum	#
#
if(is.element(design, c("triad", "cc.triad"))){
	.mm.m <- .m1.dum + .m2.dum #
	if(xchrom) .mm.f <- .f2.dum
	else .mm.f <- .f1.dum + .f2.dum #
	.mm.parents <- .mm.m + .mm.f # FOR ESTIMATING ALLELE FREQUENCIES
}
#
if(design == "cc"){
	.mm.parents <- .mm.c	# FOR ESTIMATING ALLELE FREQUENCIES
}
#
## LET THE FORMULA APPLY ONLY FOR CASES: (WARNING! RARE DISEASE ASSUMPTION HERE!)
if(design == "cc.triad"){
	#cat("Rare disease assumption made...\n")
	.mm.m <- .mm.m * .mf$cc
	.mm.f <- .mm.f * .mf$cc
	.mm.c <- .mm.c * .mf$cc
}
if(design == "cc"){
	.mm.c <- .mm.c * .mf$cc
}
if((design == "triad") & xchrom){
	.mm.girls.x.dum <- .f2.dum * .mf$sex # EXTRA EFFECT FOR GIRLS
	.mm.c <- .m2.dum + .mm.girls.x.dum # ASSUME A DOSE-RESPONSE FOR BOTH SEXES!
}
#
## SET APPROPRIATE COLUMN NAMES:
if(is.element(design, c("triad", "cc.triad"))){
	dimnames(.mm.m)[[2]] <- paste("m", 1:.A, sep = "")
	dimnames(.mm.f)[[2]] <- paste("f", 1:.A, sep = "")
}
dimnames(.mm.c)[[2]] <- paste("c", 1:.A, sep = "")
dimnames(.mm.parents)[[2]] <- paste("mf", 1:.A, sep = "")
#
##
if(F){
## DETTE KAN JO KANSKJE KOMME INN SOM DOBBELTDOSE FOR "cc",
## HVIS ANT. HAPLOTYPER ER STØRRE ENN 2 (ELLER HVA SOM TRENGS?)
	.cc <- NULL # BARE TULL, MEN FOR AA UNNGAA FEILMELDING I R CMD check
	.mm.c.dobbeldose <- (.mm.c == 2) + 0	# MODEL: R^2*Rstar
	.mm.c.dobbeldose <- .mm.c.dobbeldose * .cc$cc
	dimnames(.mm.c.dobbeldose)[[2]] <- paste("cdd", 1:.A, sep = "")	
}

#
## THE FOLLOWING IS ONLY NECESSARY IF REPONSE IS FREE
if(response == "free"){
	## SPECIAL CODING FOR HOMOZYGOTES:
	.mm.cdd <- (.mm.c == 2) + 0	# MODEL: R^2*Rstar
	if(design == "cc.triad"){
		.mm.cdd <- .mm.cdd * .mf$cc
	}
	# SET APPROPRIATE COLUMN NAMES:
	dimnames(.mm.cdd)[[2]] <- paste("cdd", 1:.A, sep = "")	#
	#
	if(maternal){
		.mm.mdd <- (.mm.m == 2) + 0	# MODEL: R^2*Rstar FOR MAT. EFF
		if(design == "cc.triad"){
			.mm.mdd <- .mm.mdd * .mf$cc			
		}
		# SET APPROPRIATE COLUMN NAMES:
		dimnames(.mm.mdd)[[2]] <- paste("mdd", 1:.A, sep = "")	
	}# END if(maternal)
}# END if(response == "free")













# SET UP FINAL DESIGN MATRIX, REMOVING REDUNDANT ref.cat-ROW.
# NOTE: DIALLELIC SITUATION REQUIRES THE REMOVAL OF ONE DOUBLE DOSE PARAMETER:

###.design.matrix <- cbind(.o) # cbind GETS THE NAMES RIGHT
if(is.element(design, c("cc.triad", "cc"))){
	.design.matrix <- .mf[,"cc", drop = F]
}
if(design == "triad"){
	if(xchrom) .design.matrix <- .mf[,"sex", drop = F]
	else .design.matrix <- NULL
}
.design.matrix <- cbind(.design.matrix, .mm.parents)

if(response != "simple"){
	.mm.c <- .mm.c[,  - ref.cat, drop = F]
	.design.matrix <- cbind(.design.matrix, .mm.c)
	if(response == "free"){
		if(.A == 2) .mm.cdd <- .mm.cdd[,  - ref.cat, drop = F]
		.design.matrix <- cbind(.design.matrix, .mm.cdd)
	}
	if(maternal){
		.mm.m <- .mm.m[,  - ref.cat, drop = F]
		.design.matrix <- cbind(.design.matrix, .mm.m)
		if(response == "free"){
			if(.A == 2) .mm.mdd <- .mm.mdd[,  - ref.cat, drop = F]
			.design.matrix <- cbind(.design.matrix, .mm.mdd)
		}
	}
}

.design.matrix <- as.data.frame(.design.matrix)

if(test & (response == "mult")){

	# DETTE ER EN EGEN TEST FOR AA PROEVE AA MATCHE RESULTATET FRA XLRT,
	# VED EN KJOERING
	# jj <- haplin(testdX1, xchrom = T, sex = 1, n.vars = 2, markers = 1, verbose = T, printout = F, response = "mult", reference = 2, use.missing = T)
	# OG DERETTER SE PAA  exp(jj$result$result$coef)
	# MERK AT f.tri.glm KALLER f.make.design MED test = T

	.sex <- .design.matrix$sex
	.c1 <- .design.matrix$c1
	#
	.a1 <- (.c1 == 1) & (.sex == 0)
	.a2 <- (.c1 == 0) & (.sex == 1)
	.a3 <- (.c1 == 1) & (.sex == 1)
	.a4 <- (.c1 == 2) & (.sex == 1)
	#
	.design.matrix$a1 <- .a1 + 0
	.design.matrix$a2 <- .a2 + 0
	.design.matrix$a3 <- .a3 + 0
	.design.matrix$a4 <- .a4 + 0
	#
	.design.matrix$sex <- .design.matrix$c1 <- NULL
	#
	.msum <- .m1.dum[,2] + .m2.dum[,2]
	.fsum <- .f2.dum[,2]
	.mu <- paste(.msum, .fsum, sep = "-")
	.mu <- model.matrix(~ -1 + factor(.mu))
	colnames(.mu) <- paste("mu", seq(length.out = dim(.mu)[2]), sep = "")
	.design.matrix$mf1 <- .design.matrix$mf2 <- NULL
	.design.matrix <- cbind(.mu, .design.matrix)

	### tull <<- .design.matrix
	#.design.matrix <- cbind(.m1.dum, .m2.dum, .f2.dum, .mu, .design.matrix)
	#stop()
}



#.formula <- formula(paste(c(".o ~ -1 ", names(.design.matrix)), collapse = "+"))
#print(.formula)

return(.design.matrix)

}
