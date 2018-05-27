print.haplin.data <- function( x, ... ){
	message( "This is raw genetic data read in through genDataRead.\n
It contains the following parts:\n  ", paste( names( x ), collapse = ", " ), "\n
with following dimensions:" )
	if( !is.null( x$cov.data ) ){
		message( "
  - covariate variables = ", paste( colnames( x$cov.data ), collapse = ", " ), "
      (total ", ncol( x$cov.data ), " covariate variables),
" )
	} else {
		message( "  - no covariate variables," )
	}
	
	format <- x$aux$info$filespecs$format
	if( format == "ped" ){
		fam.size <- 2
	} else { #format == "haplin"
		fam.size <- 6
	}

	message( "  - number of markers = ", sum( unlist( lapply( x$gen.data, ncol ) ) )/fam.size, "," )
	message( "  - number of data lines = ", nrow( x$gen.data[[ 1 ]] ) )
}
