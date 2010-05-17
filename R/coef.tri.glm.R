coef.tri.glm <- function(object, ...){
##
## EXTRACTS COEFFICIENTS AND THE VARIANCE-COVARIANCE MATRIX FROM A tri.glm OBJECT
##
#
## EXTRACT glm PART
.res <- object$result
#
.summary <- summary.glm(.res, correlation = F)
.coef <- .summary$coefficients[, 1] # WORKS FOR BOTH SPLUS AND R
.cov <- .summary$cov.unscaled	##
#
## REPLACE "RAW" GLM COVARIANCE WITH EM-CORRECTED, IF EXISTS:
.cov.correct <- attr(object, "cov.correct")
if(!is.null(.cov.correct)) .cov <- .cov.correct
#
## USE RESAMPLED MATRIX IN PREFERENCE TO THE OTHERS
.cov.resamp <- attr(object, "cov.resamp")
if(!is.null(.cov.resamp)) .cov <- .cov.resamp
#
##
.ut <- list(coef = .coef, cov = .cov)
#
return(.ut)

}
