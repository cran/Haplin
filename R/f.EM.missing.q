"f.EM.missing" <- function(data, maternal, ref.cat, design = "triad", xchrom, response = "free", max.EM.iter, x = F, verbose = T, suppressEmWarnings = F)
{
## PERFORMS THE ESSENTIALS OF THE EM-ALGORITHM, INCLUDING HANDLING AMBIGUOUS HAPLOTYPES
## AND MISSING GENOTYPE DATA
#


#### WARNINGS: ##############################
# MERK DENNE!!
## FOR GLM: abs(new.deviance-old.deviance)/(old.deviance + epsilon) < epsilon
## HVOR epsilon = 0.0001
	f.vis("Sjekk kriteriene!", vis = F)

#
#### PREPARE: ############################
	.n.haplo <- sum(attr(data, "selected.haplotypes"))
###	.tmpfreq <- rep(1, length(.agg.freq))
if((design == "triad") & !xchrom){
	.tmpfreq <- rep(1, .n.haplo^4) # ALL POSSIBLE GENOTYPES FOR TRIAD
}
if((design == "triad") & xchrom){
	.tmpfreq <- rep(1, 2 * .n.haplo^3) # ALL POSSIBLE GENOTYPES FOR TRIAD, BOTH MALES AND FEMALES
}
if(design == "cc"){
	if(xchrom)stop("Not implemented!")
	.tmpfreq <- rep(1, 2 * .n.haplo^2) # ALL POSSIBLE GENOTYPES FOR BOTH CASES AND CONTROLS
}
if(design == "cc.triad"){
	if(xchrom)stop("Not implemented!")
	.tmpfreq <- rep(1, 2 * .n.haplo^4) # ALL POSSIBLE GENOTYPES FOR TRIAD, BOTH CASES AND CONTROLS
}


###	.agg.freq <- agg.data$freq	#
###	.amb <- data$amb
###	.freqsum <- f.groupsum(X = .freq, INDICES = .amb)
###	.denominator <- f.groupsum(X = rep(1, length(.amb)), INDICES = .amb)	# COULD USE tmp FROM ORIGINAL MATRIX, BUT f.thin HAS BEEN USED IN THE MEANTIME...
# COMPUTE STARTING FREQUENCIES FOR THE EM BY TAKING AVERAGE OVER AMBIGUOUS CATEGORIES:
###	.tmpfreq <- .freqsum/.denominator	#
##	.tmpfreq <- .freq
#
#### INITIALIZE EM LOOP: ##################
.deviance <- numeric(max.EM.iter)	#
i <- 1	#
.EM.conv <- F
#### EM LOOP: ###########################
repeat{
# M-STEP:
	f.vis(sum(.tmpfreq), vis = F)

	.res <- f.tri.glm(.tmpfreq, maternal = maternal, ref.cat = ref.cat, response = response, design = design, xchrom = xchrom)	#

	.deviance[[i]] <- .res$result$deviance	# THIS IS TWICE THE MINUS LOG-LIKELIHOOD OF THE GLM
	.coef.new <- .res$result$coefficients
	#
	## PRINT INTERMEDIATE RESULTS (IF REQUESTED):
	if(verbose) {
		#	cat("-----------------------------------------------\n")
			cat("EM iter:", sprintf("%-4.0f", round(i)), "|")	#
		#	cat("-----------------------------------------------\n")
		#	cat("GLM deviance:", format(.deviance[[i]], width = 6, scientific = F, digits = 6), "|")
			cat("GLM deviance:", sprintf("%-12.6g", .deviance[[i]]), "|")
		#	cat("Coefficients:", format(.coef.new, width = 6, scientific = F, digits = 6), "\n")
			cat("Coefficients:", sprintf("%-12.6g", .coef.new), "\n")
		#	cat("-----------------------------------------------\n")
	}# END if(verbose)
	else{
		# cat(".") # SKRUDD AV, MIDLERTIDIG
	}
	#
	## STOP WHEN CONVERGED
	if(i > 1) {
		#
		## STOPPING CRITERIA:
		.crit.deviance <- abs(.deviance[i] - .deviance[i - 1]) < 2e-006
		.crit.coef <- max(abs(.coef.new - .coef.old)) < 1e-006
		if(.crit.deviance & .crit.coef){
			.EM.conv <- T
			break
			}
	}# END if(i > 1)
	#
	## BREAK OFF, WITH WARNING, IF max.EM.iter IS EXCEEDED:
	if(i >= max.EM.iter) {
		if(!suppressEmWarnings) warning("Maximum number of EM iterations reached!\n Convergence not yet obtained. Setting max.EM.iter higher may help.")
		break
	}
#
## UPDATING VALUES, E-STEP:  
##
		i <- i + 1
		.coef.old <- .coef.new
		.pred <- .res$pred	#
		.tmpfreq <- f.redistribute(pred = .pred, data = data, design = design, xchrom = xchrom)
		

###		.predsum <- f.groupsum(X = .pred, INDICES = .amb)	#
##  .tmpfreq <- .pred/.predsum * .freqsum
###		.tmpfreq <- ifelse(.predsum > 0, .pred/.predsum * .freqsum, 0)

#
## CHECK CONSISTENCY OF NEW VALUES:
		if(any(is.na(.pred))) {
			warning("Missing in predicted values.... There may be a problem with the EM algorithm!")
			.tmpfreq[is.na(.tmpfreq)] <- 1e-005	# FORSIKTIG HER!!!
		}
		if(sum(is.na(.tmpfreq)) > 0) {
			stop("Missing values in EM-updated frequencies!")
		}
next
}# END REPEAT
#
if(x){
	## ADD DESIGN MATRIX IN LAST STEP, IF REQUESTED:
	.res <- f.tri.glm(.tmpfreq, maternal = maternal, ref.cat = ref.cat, response = response, design = design, xchrom = xchrom, x = x)#
} # WARNING: REQUIRES .tmpfreq TO BE UNCHANGED AFTER LAST COMPUTATION OF .res!
#
attr(.res, "iter.used") <- i
attr(.res, "EM.conv") <- .EM.conv
#tull <<- .res
return(.res)
}
