f.expand.alleles <- function( data, max.allele, nmarkers, nalleles ){
	for(i in max.allele:1){# FOR EACH OF THE X ALLELES
		for(j in 1:nmarkers){
			# HOW MANY LINES TO EXPAND EACH MISSING INTO:
			.nexpand <- rep(1, dim(data)[1])
			.nexpand[is.na(data[,i]) & (data[,"ind.marker"] == j)] <- nalleles[j]
			# CREATE THE INDEX FOR EXPANSION AND EXPAND:
			.ind.expand <- rep(seq(along = .nexpand), .nexpand)
			data <- data[.ind.expand, , drop = F]
			# REPLACE THE NAs IN THE EXPANDED DATA WITH SEQUENCE OF ALL POSSIBLE ALLELES AT MARKER:
			if( max.allele == 3 & i == 3 ){
				data[is.na(data[,3]) & (data[,"ind.marker"] == j),3:4] <- 1:nalleles[j]
			} else {
				data[is.na(data[,i]) & (data[,"ind.marker"] == j),i] <- 1:nalleles[j]
			}
		}
	}
	return( data )
}
