print.haplin.ready <- function( x, ... ){
	message( "
		This is preprocessed data, ready for haplin analysis.\n
		It contains the following parts:\n  ", paste( names( x ), collapse = ", " ), "\n
		with following dimensions:
	" )
	if( !is.null( x$cov.data ) ){
		message( "  - number of covariate variables = ", ncol( x$cov.data ), "," )
	} else {
		message( "  - no covariate variables," )
	}
	message( "  - number of markers = ", length( x$aux$marker.names ), "," )
	message( "  - number of individuals/families = ", nrow( x$gen.data[[ 1 ]] ) )
}
