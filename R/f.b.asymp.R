f.b.asymp <- function(data, p, ncells, norig, orig, info){
	## Set up ambiguity matrix
	## Then multiply with p and standardize
	.a <- f.a.asymp(data, ncells = ncells, norig = norig, orig = orig, info = info)
	.ap <- .a*p
	.b <- t(t(.ap)/colSums(.ap))
	return(.b)
}
