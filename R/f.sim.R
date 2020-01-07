f.sim <- function(.prob, size, nall, .nloci, xchrom, .grid){
##
## SAMPLE FROM MULTINOMIAL DISTRIBUTION
##	.sim.rownum <- rMultinom(matrix(.prob, nrow = 1), m = n.cases)
	.sim.rownum <- suppressWarnings(sample(length(.prob), size = size, replace = T, prob = .prob))
	#
	## FIND ALLELES CORRESPONDING TO GRID ROW NUMBERS
	if(xchrom){
		.nrows <- dim(.grid)[1]/2
		.girls <- (.sim.rownum > .nrows)
		.tmp.boys <- .sim.rownum[!.girls]
		.tmp.girls <- .sim.rownum[.girls]
		.tmp.girls <- .tmp.girls - .nrows
		if(length(.tmp.boys!=0)) .alleles.boys <- f.pos.to.haplocomb(A = nall, pos = .tmp.boys, fam = "mfx")
		if(length(.tmp.girls!=0)) .alleles.girls <- f.pos.to.haplocomb(A = nall, pos = .tmp.girls, fam = "mfx")
		###.sex1 <- .sex[.sim.rownum]
		.sex1 <- c(rep(1, length(.tmp.boys)), rep(2, length(.tmp.girls)))
		###.alleles <- rbind(.alleles.boys, .alleles.girls)
	} else{
		.alleles <- f.pos.to.haplocomb(A = nall, pos = .sim.rownum)
	}
	#
	## ADD COLUMNS WITH CHILD GENOTYPES
	if(xchrom){
		if(length(.tmp.boys!=0)) .names <- dimnames(.alleles.boys)[[2]]
		else .names <- dimnames(.alleles.girls)[[2]]

		.names <- matrix(.names, nrow = .nloci)
		
		if(length(.tmp.boys!=0)){
			.ind.boys <- c(1,2,3,3,2,2)
			.names.boys <- .names[, .ind.boys]
			.names.boys <- as.vector(t(.names.boys))
			.alleles.boys <- .alleles.boys[,.names.boys]
			.alleles <- .alleles.boys
		}
		
		if(length(.tmp.girls!=0)){
		.ind.girls <- c(1,2,3,3,2,3)		
			.names.girls <- .names[, .ind.girls]
			.names.girls <- as.vector(t(.names.girls))
			.alleles.girls <- .alleles.girls[,.names.girls]
			.alleles <- .alleles.girls
		}
		if(length(.tmp.boys!=0) & length(.tmp.girls!=0)).alleles <- rbind(.alleles.boys, .alleles.girls)
	}else{
		.names <- dimnames(.alleles)[[2]]
		.names <- matrix(.names, nrow = .nloci)
		.ind <- c(1,2,3,4,2,4)
		.names <- .names[, .ind]
		.names <- as.vector(t(.names))
		#
		.alleles <- .alleles[,.names]
	}
	#
	if(is.vector(.alleles)) .alleles <- as.matrix(t(.alleles))
	## RANDOMIZE SEQUENCE OF ALLELES AND PREPARE FOR WRITING TO DISK
	.all.paste <- vector(dim(.alleles)[2]/2, mode = "list")
	for(j in seq(along = .all.paste)){
		.all.paste[[j]] <- f.rand.geno(.alleles[,2*j-1], .alleles[,2*j])
	}

	.all.paste <- as.data.frame(.all.paste)
	names(.all.paste) <- seq(along = .all.paste)

	if(xchrom) .ut <- cbind(sex = .sex1, .all.paste)
	else .ut <- .all.paste
	#
	return(.ut)	
}	
