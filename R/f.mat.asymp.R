f.mat.asymp <- function(b, w, ncells){
	## Sum up diag(b) - b bT
	## Apply weights w[i]
	## Can probably be sped up, cf f.var.covar
	## WARNING!: w should match .orig!!
	.norig <- ncol(b) # equals to .ncells, but may expand later
	if(missing(w)) w <- rep(1, .norig)
	.matsum <- matrix(0, nrow = ncells, ncol = .norig)
	for(i in 1:.norig){
		.b <- b[,i]
		.mat <- diag(.b) - .b %*% t(.b)
		.matsum <- .matsum + w[i]*.mat
	}
	return(.matsum)
}
