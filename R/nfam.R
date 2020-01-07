#' Count the number of families in the data
#'
#' This is a help function to count the number of families in an object read in with
#'  \link{genDataRead} (or loaded with \link{genDataLoad}). Note: it is assumed that
#'  the study design is either 'triad' or 'cc.triad'.
#'
#' @param data.in The data read in by \link{genDataRead}.
#'
#' @return How many families (integer).
#'
nfam <- function( data.in ){
	# check if input data is in correct format
	if( !is( data.in, "haplin.data" ) ||
	  !all( names( data.in ) == .haplinEnv$.haplin.data.names ) ){
		stop( "The input data is not in the correct format!", call. = FALSE )
	}

	format <- data.in$aux$info$filespecs$format
	fam.id.colname <- get( ".cov.data.colnames", envir = .haplinEnv )[1]
	# if it's a PED data or haplin with "fam.id" column in the covariate data
	if( !is.null( data.in$cov.data ) &
		( fam.id.colname %in% colnames( data.in$cov.data ) ) ){
		all.families <- length( unique( data.in$cov.data[, fam.id.colname ] ) )
	} else if( format == "haplin" ){
		all.families <- nrow( data.in$gen.data[[1]] )
	} else {
		stop( "Cannot make out the number of families. Are you sure your data is in the correct format?", .call = FALSE )
	}

	return( all.families )
}
