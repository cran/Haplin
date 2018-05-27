f.extract.freq <- function( data, alleles.orig, design ){
	freq.x <- NULL
	nas.x <- NULL
	if( !is.null( data$cov.data ) ){
		## GET FREQUENCY COUNT FOR COVARIATE DATA
		tmp <- lapply( 1:ncol( data$cov.data ), function(i) f.freq.table( unlist( data$cov.data[,i] ) ) )
		names( tmp ) <- colnames( data$cov.data )
		freq.x <- lapply( tmp, function(x) x$tab )
		nas.x <- sapply( tmp, function(x) x$nas )
	}
	
	if( design %in% c( "triad", "cc.triad" ) ) {
		tt <- 6
	}else if( design == "cc" ) tt <- 2
	marker.seq <- 1:( ncol( data$gen.data ) / tt )
	
	alleles.freq <- sapply( marker.seq, function( i ){
		list( f.freq.table( data$gen.data[ ,( (i-1)*tt + 1 ):( i*tt ) ] ) )
	} )
	names( alleles.freq ) <- names( alleles.orig )
	alleles.new <- lapply( 1:length( alleles.freq ), function(x){
		freq.tmp <- alleles.freq[[x]]$tab
		names( freq.tmp ) <- names( alleles.orig[[x]] )[ as.numeric( names( freq.tmp ) ) ]
		return( freq.tmp )
	})
	names( alleles.new ) <- names( alleles.orig )
	alleles.nas.new <- lapply( alleles.freq, function(x) x$nas )

	return( list( variables = freq.x, variables.nas = nas.x, alleles = alleles.new, alleles.nas = alleles.nas.new ) )
}
