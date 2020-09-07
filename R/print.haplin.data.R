print.haplin.data <- function( x, ... ){
	cat( "This is raw genetic data read in through genDataRead.\n
It contains the following parts:\n  ", paste( names( x ), collapse = ", " ), "\n
with following dimensions:" )
	if( !is.null( x$cov.data ) ){
		cat( "
  - covariate variables = ", paste( colnames( x$cov.data ), collapse = ", " ), "
      (total ", ncol( x$cov.data ), " covariate variables),
" )
	} else {
		cat( "  - no covariate variables,\n" )
	}
	
	format <- x$aux$info$filespecs$format
	if( format == "ped" ){
		fam.size <- 2
	} else { #format == "haplin"
		fam.size <- 6
	}

	cat( "  - number of markers = ", sum( unlist( lapply( x$gen.data, ncol ) ) )/fam.size, ",\n" )
	cat( "  - number of data lines = ", nrow( x$gen.data[[ 1 ]] ), "\n" )
}
