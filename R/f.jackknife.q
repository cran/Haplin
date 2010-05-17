f.jackknife.new <- function(data, maternal, ref.cat, verbose = F, use.EM, max.EM.iter, ...){

### TIL NY VARIANT ### f.jackknife <- function(data, design.matrix, verbose = F...){

## USES JACKKNIFE RESAMPLING TO COMPUTE STANDARD ERRORS AND/OR CONFIDENCE INTERVALS
## USES A PRECOMPUTED DESIGN MATRIX
##
#
## PREPARE:
.ind <- unique(data$ind)
.coef.jack <- vector(length(.ind), mode = "list")
#
## JACKKNIFE LOOP:
for (i in seq(along = .ind)){
	.data.tmp <- data[-which(data$ind == .ind[i]),] # JACKKNIFE, ONE FOR EACH TRIAD
	cat("\nJackknife sample no.", i, "of", length(.ind), "total\n")
	
###		cat("kan unngaas! Kan subtrahere fra data.agg\n")
###		.data.agg.tmp <- f.sum.and.expand(.data.tmp)		


### TIL NY VARIANT ###	.coef.jack[[i]] <- f.EM.test(data = .data.tmp, design.matrix = design.matrix, verbose = verbose, ...)$result$coefficients
	if(use.EM){
	.coef.jack[[i]] <- f.EM.missing(.data.tmp, maternal = maternal, xchrom = F, max.EM.iter = max.EM.iter, verbose = verbose, ref.cat = ref.cat)$result$coefficients}
	else{
		stop("Sjekk denne! Her skal det vel vaere .data.agg.tmp?")
		.coef.jack[[i]] <- f.tri.glm(.data.tmp$freq, maternal = maternal, ref.cat = ref.cat)$result$coefficients}
}# END for i
#
## SET IN MATRIX FORM:
.check <- sapply(.coef.jack, length)
if(any(.check != .check[1])) stop("Problem: differing number of coefficients in jackknife!\n")
.coef.jack.mat <- matrix(unlist(.coef.jack), nrow = length(.coef.jack[[1]]), dimnames = list(names(.coef.jack[[1]]), NULL)) #
#
## COMPUTE WEIGHTED MEAN AND VAR-COVAR MATRIX:
###	.posfreq <- data$freq[.ind]
.posfreq <- rep(1, length(.ind)) # MODIFISERT MIDLERTIDIG

.mean <- as.numeric(.coef.jack.mat %*% .posfreq / sum(.posfreq))
.cent.coef <- .coef.jack.mat - .mean
.cov <- .cent.coef %*% (.posfreq * t(.cent.coef)) #
#
#
return(list(mean = .mean, cov = .cov))









if(F){

	##### KEEP THIS, FOR NOW. WILL INCLUDE FREQUENCIES LATER #####

	## PREPARE:
		.ind <- which(data$freq >= 1)
		.coef.jack <- vector(length(.ind), mode = "list") #
	#
	## JACKKNIFE LOOP:
		for (i in seq(along = .ind)){
			.data.tmp <- data
			.data.tmp[.ind[i],"freq"] <- .data.tmp[.ind[i],"freq"] - 1 # JACKKNIFE, ONE FOR EACH CELL
			cat("\nJackknife sample no.", i, "of", length(.ind), "total\n")
	### TIL NY VARIANT ###	.coef.jack[[i]] <- f.EM.test(data = .data.tmp, design.matrix = design.matrix, verbose = verbose, ...)$result$coefficients
			if(use.EM){
			.coef.jack[[i]] <- f.EM.missing(.data.tmp, maternal = maternal, max.EM.iter = max.EM.iter, verbose = verbose, ref.cat = ref.cat)$result$coefficients}
			else{.coef.jack[[i]] <- f.tri.glm(.data.tmp$freq, maternal = maternal, ref.cat = ref.cat)$result$coefficients}
		}# END for i
	#
	## SET IN MATRIX FORM:
		.check <- sapply(.coef.jack, length)
		if(any(.check != .check[1])) stop("Problem: differing number of coefficients in jackknife!\n")
		.coef.jack.mat <- matrix(unlist(.coef.jack), nrow = length(.coef.jack[[1]]), dimnames = list(names(.coef.jack[[1]]), NULL)) #
	#
	## COMPUTE WEIGHTED MEAN AND VAR-COVAR MATRIX:
		.posfreq <- data$freq[.ind]
		.mean <- as.numeric(.coef.jack.mat %*% .posfreq / sum(.posfreq))
		.cent.coef <- .coef.jack.mat - .mean
		.cov <- .cent.coef %*% (.posfreq * t(.cent.coef)) #
	#
	#
		return(list(mean = .mean, cov = .cov))
}

}
