coef.haplin <- function(object, ...){
##
## EXTRACTS COEFFICIENTS AND THE VARIANCE-COVARIANCE MATRIX FROM A haplin OBJECT
##
.res <- object$result

.ut <- coef.tri.glm(.res)

return(.ut)


}
