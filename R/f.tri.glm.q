"f.tri.glm" <- function(observed, maternal = F, ref.cat, response = "free", design = "triad", ...)
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
# OBSERVED FREQUENCY, ARRANGED ACCORDING TO GRID:
	.o <- observed	#
if(design == "triad"){
	.A <- round(length(.o)^(1/4))	# NUMBER OF ALLELES. NOTE: ASSUMES DATA FROM FULL GRID!
	if(abs(.A^4 - length(.o)) > 1e-007)
		stop("Data vector of incorrect length! Assumes data from\n  full grid.\n")
}
if(design == "cc.triad"){
	.A <- round((length(.o)/2)^(1/4))	# NUMBER OF ALLELES. NOTE: ASSUMES DATA FROM FULL GRID!
	if(abs(2 * .A^4 - length(.o)) > 1e-007)
		stop("Data vector of incorrect length! Assumes data from\n  full grid.\n")
}
if(design == "cc"){
	.A <- round((length(.o)/2)^(1/2))	# NUMBER OF ALLELES. NOTE: ASSUMES DATA FROM FULL GRID!
	if(abs(2 * .A^2 - length(.o)) > 1e-007)
	stop("Data vector of incorrect length! Assumes data from\n  full grid.\n")
}
	.ntri <- sum(.o)	# NUMBER OF TRIADS (OR INDIVIDUALS, FOR CASE-CONTROL)
#
# CREATE GRID OF PARENTAL (OR FETAL) GENOTYPES (AND/OR POSSIBLY CASE-CONTROL STATUS)
.tmp <- factor(1:.A)
if(design == "triad"){
	.mf <- expand.grid(m1 = .tmp, m2 = .tmp, f1 = .tmp, f2 = .tmp)	#
}
if(design == "cc.triad"){	
	.mf <- expand.grid(m1 = .tmp, m2 = .tmp, f1 = .tmp, f2 = .tmp, cc = c(0,1)) ## HERE, 0 AND 1 (CONTROL AND CASE, RESP.) ARE USED AS DUMMY. ORIGINAL CODING (OUTSIDE f.tri.glm) SHOULD BE 1 (CONTROL) AND 2 (CASE). IN INPUT FILE, THE LARGEST IS ALWAYS THE CASE
}
if(design == "cc"){	
	.mf <- expand.grid(m2 = .tmp, f2 = .tmp, cc = c(0,1)) ## HERE, 0 AND 1 (CONTROL AND CASE, RESP.) ARE USED AS DUMMY. ORIGINAL CODING (OUTSIDE f.tri.glm) SHOULD BE 1 (CONTROL) AND 2 (CASE)
}
#
## CREATE "DUMMY" (COUNTING) VARIABLES FOR PATERNAL GENOTYPES AND FOR CHILDREN:
.mm.c <- model.matrix( ~ -1 + m2, data = .mf) + model.matrix( ~ -1 + f2, data = .mf)	#
#
if(is.element(design, c("triad", "cc.triad"))){
	.mm.m <- model.matrix( ~ -1 + m1, data = .mf) + model.matrix( ~ -1 + m2, data = .mf)	#
	.mm.f <- model.matrix( ~ -1 + f1, data = .mf) + model.matrix( ~ -1 + f2, data = .mf)	#
	.mm.parents <- .mm.m + .mm.f	# FOR ESTIMATING ALLELE FREQUENCIES
}
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

.design.matrix <- cbind(.o) # cbind GETS THE NAMES RIGHT
if(is.element(design, c("cc.triad", "cc"))){
	.design.matrix <- cbind(.design.matrix, .mf[,"cc", drop = F])
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

.formula <- formula(paste(c(".o ~ -1 ", names(.design.matrix)[-1]), collapse = "+"))

f.vis(apply(.design.matrix$.o * .design.matrix, 2, sum), vis = F)

#
# ESTIMATION:
f.vis("iteration controls....?", vis = F)
f.vis("Should rather use glm.fit here...", vis = F)

.res <- suppressWarnings(glm(.formula, family = poisson, data = .design.matrix, ..., maxit = 20))	#
# 
#
# FREQUENCY PREDICTION (NOTE: COULD ALSO HAVE USED FITTED VALUES IN OBJECT): 
.pred <- predict(.res, type = "response")	#
if(abs(sum(.pred) - .ntri) > 0.001 * .ntri) stop("Potential problem in prediction!")	#
#
#
# ADDING INFORMATION TO OUTPUT:
.out <- list(result = .res, pred = .pred, nall = .A, ntri = .ntri, ref.cat = ref.cat, maternal = maternal, design = design, orig.call = sys.call(), date = date())
class(.out) <- "tri.glm"	#
return(invisible(.out))	# 
}
