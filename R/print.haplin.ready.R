print.haplin.ready <- function( x, ... ){
	cat( "
		This is preprocessed data, ready for haplin analysis.\n
		It contains the following parts:\n  ", paste( names( x ), collapse = ", " ), "\n
		with following dimensions:
	" )
	if( !is.null( x$cov.data ) ){
		cat( "  - number of covariate variables = ", ncol( x$cov.data ), ",\n" )
	} else {
		cat( "  - no covariate variables,\n" )
	}
	cat( "  - number of markers = ", length( x$aux$marker.names ), ",\n" )
	cat( "  - number of individuals/families = ", nrow( x$gen.data[[ 1 ]] ), "\n" )
}
