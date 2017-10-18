f.a.asymp <- function(data, ncells, norig, orig, info){
	## Set up ambiguity matrix
	##
	.a <- matrix(0, nrow = ncells, ncol = norig)
	dimnames(.a) <- list(gridln = 1:ncells, origln = orig)
	# Loop over original lines
	for(i in seq(along = orig)){
		.tmpd <- data[data$orig.lines == orig[i], , drop = F]
		.pos <- f.pos.match(data = .tmpd, info = info)
		.a[.pos,i] <- 1
	}
	return(.a)
}