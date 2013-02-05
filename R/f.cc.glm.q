"f.cc.glm"<-
function(observed, ref.cat, response = "free", ...)
{
#
# THE PROGRAM ESTIMATES EFFECTS OF SEVERAL ALLELES IN A CASE-CONTROL DESIGN
# ASSUMING HARDY-WEINBERG EQUILIBRIUM AND RARE DISEASE
# IMPORTANT: OBSERVED FREQUENCIES (observed) MUST COME FROM A COMPLETE
# GRID IN APPROPRIATE ORDER! THE ORDER IS BY FIRST ALLELE, SECOND ALLELE
# AND FINALLY CONTROL (1) AND CASE (2), SO THAT TOP HALF OF THE DATA ARE
# CONTROLS
# THE ... PASSES ON OTHER ARGUMENTS, LIKE start, TO THE GLM
#
# OBSERVED FREQUENCY, ARRANGED ACCORDING TO GRID:
	.o <- observed	#
	.A <- round((length(.o)/2)^(1/2))	# NUMBER OF ALLELES. NOTE: ASSUMES DATA FROM FULL GRID!
	if(abs(2 * .A^2 - length(.o)) > 1e-007)
		stop("Data vector of incorrect length! Assumes data from\n  full grid.\n")
	.nind <- sum(.o)	# NUMBER OF INDIVIDUALS
#
# CREATE GRID OF GENOTYPES AND CASE-CONTROL STATUS
	.tmp <- factor(1:.A)
	.cc <- expand.grid(c1 = .tmp, c2 = .tmp, cc = c(0,1)) ## HERE, 0 AND 1 (CONTROL AND CASE, RESP.) ARE USED AS DUMMY. ORIGINAL CODING (OUTSIDE f.cc.glm) SHOULD BE 1 (CONTROL) AND 2 (CASE)

















#
# CREATE "DUMMY" (COUNTING) VARIABLES FOR GENOTYPES:
	.mm.c <- model.matrix( ~ -1 + c1, data = .cc) + model.matrix( ~ -1 + c2, data = .cc)	#
	.mm.p <- .mm.c	# FOR ESTIMATING ALLELE FREQUENCIES




















#
# SPECIAL CODING FOR HOMOZYGOTES: 
	.mm.c.dobbeldose <- (.mm.c == 2) + 0	# MODEL: R^2*Rstar
# NOTE: THE FOLLOWING CODINGS FOR THE DOUBLE DOSE WILL DEPEND ON REFERENCE CATEGORY!:
## .mm.c <- .mm.c - 2 * .mm.c.dobbeldose # MODEL: 1*Rstar
## .mm.c <- .mm.c - 1 * .mm.c.dobbeldose # MODEL: R*Rstar
#
## LET THE FORMULA APPLY ONLY FOR CASES: (WARNING! RARE DISEASE ASSUMPTION HERE!)
cat("Rare disease assumption made...\n")
.mm.c <- .mm.c * .cc$cc
.mm.c.dobbeldose <- .mm.c.dobbeldose * .cc$cc
#
# SET APPROPRIATE COLUMN NAMES:
	dimnames(.mm.p)[[2]] <- paste("p", 1:.A, sep = "")
	dimnames(.mm.c)[[2]] <- paste("c", 1:.A, sep = "")
	dimnames(.mm.c.dobbeldose)[[2]] <- paste("cdd", 1:.A, sep = "")	#
#
# SET UP FINAL DESIGN MATRIX, REMOVING REDUNDANT ref.cat-ROW.
# NOTE: DIALLELIC SITUATION REQUIRES THE REMOVAL OF ONE DOUBLE DOSE PARAMETER:

	if(.A == 2) .design.matrix <- as.data.frame(cbind(.o, .cc[,"cc", drop = F], .mm.p, .mm.c[,  - ref.cat, drop = F], .mm.c.dobbeldose[,  - ref.cat, drop = F])) 
	else .design.matrix <- as.data.frame(cbind(.o, .cc[,"cc", drop = F], .mm.p, .mm.c[, -ref.cat, drop = F], .mm.c.dobbeldose)) # cbind HELPS PRODUCING THE RIGHT NAMES
#




	.formula.simple <- formula(paste(c(".o ~ -1 ", names(.design.matrix)[2:(.A + 2)]), collapse = "+"))
###	.formula.mult <- formula(paste(c(".o ~ -1 ", names(.design.matrix)[2:(2 * .A)]), collapse = "+"))	#
	.formula.mult <- formula(paste(c(".o ~ -1 ", names(.design.matrix)[2:(2 * .A + 1)]), collapse = "+"))	#
	.formula.free <- formula(paste(c(".o ~ -1 ", names(.design.matrix)[-1]), collapse = "+"))	#
f.vis(.formula.simple)
f.vis(.formula.mult)
f.vis(.formula.free)

if(response == "simple") .formula <- .formula.simple
if(response == "mult") .formula <- .formula.mult
if(response == "free") .formula <- .formula.free

#
# ESTIMATION:
		f.vis("iteration controls....?")
		f.vis("Should rather use glm.fit here...")

		.res <- suppressWarnings(glm(.formula, family = poisson, data = .design.matrix, ..., maxit = 20))	#
# 
# FREQUENCY PREDICTION (NOTE: COULD ALSO HAVE USED FITTED VALUES IN OBJECT): 
	.pred <- predict(.res, type = "response")	#
	if(abs(sum(.pred) - .nind) > 0.001 * .nind) stop("Potential problem in prediction!")	#
#
#
# ADDING INFORMATION TO OUTPUT:
	.out <- list(result = .res, pred = .pred, nall = .A, ntri = .nind, ref.cat = ref.cat, maternal = F, design = "cc", orig.call = sys.call(), date = date())

	class(.out) <- "tri.glm"	#
	return(invisible(.out))	# 
}
