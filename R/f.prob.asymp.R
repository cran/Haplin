f.prob.asymp <- function(beta, design, X, k){
#
## Design matrix and probabilities
.design <- design
.beta <- beta
.X <- X
.k <- k
.prob <- exp(.X%*%.beta)
if(.design!="triad"){
	.n <- length(.prob)
	.prob.cases <- sum(.prob[(.n/2+1):.n])
	.prob.controls <- sum(.prob[1:(.n/2)])
	.beta.cc <- log(.prob.controls/(.prob.cases*.k))
	.beta[which(names(.beta)=="cc")] <- .beta.cc
	.prob <- exp(.X%*%.beta)
	.prob <- as.vector(.prob/sum(.prob))
} else .prob <- as.vector(.prob/sum(.prob))

return(.prob)

}
