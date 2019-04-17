#' Display phenotype part of data
#'
#' This is a help function to extract phenotypic (and covariate) data from an object
#' read in with \link{genDataRead} (or loaded with \link{genDataLoad}).
#'
#' @param data.in The data read in by \link{genDataRead}).
#' @param n Number of rows to display or "all" (default: 5).
#' @param from From which row to display (optional, default: from the first).
#' @param to To which row to display (optional).
#' @param sex If the sex column is part of the phenotypic information, the user can
#'    choose based on one of the categories used in this column (optional);
#'    NB: this does not combine with the 'to' and 'from' arguments.
#'
#' @return A table with phenotypic and covariate data (if any) extracted from 'data.in'.
#'

### @param ids Individual IDs to display (optional); NB: this does not combine with the
###    'to' and 'from' arguments.

showPheno <- function( data.in, n = 5, from, to, sex ){
	# check if input data is in correct format
	if( ( class( data.in ) != "haplin.data" ) || !all( names( data.in ) == .haplinEnv$.haplin.data.names ) ){
		stop( "The input data is not in the correct format!", call. = FALSE )
	}

	# check if there is any phenotypic info
	if( is.null( data.in$cov.data ) ){
		stop( "There is no phenotypic/covariate information in your data!", call. = FALSE )
	}
	
	cov.data <- data.in$cov.data
	n.given <- FALSE
	from.given <- FALSE

	if( !missing( sex ) ){
		if( !( "sex" %in% colnames( cov.data ) ) ){
			stop( "There is no 'sex' column in the dataset!" )
		}
		if( !( sex %in% as.numeric( unique( cov.data[ ,"sex" ] ) ) ) ){
			stop( paste( "Wrong sex chosen - 'sex' column has the following categories:", unique( cov.data[ ,"sex" ] ), collapse = "," ), call. = FALSE )
		}
		
		which.rows.show <- as.numeric( cov.data[ ,"sex" ] ) == sex
		data.out <- cov.data[ which.rows.show, ]
	} else {
		if( !missing( n ) ){
			n.given <- TRUE
			if( is.numeric( n ) ){
				if( n < 1 | n > nrow( data.in$cov.data ) ){
					stop( "Your specification of 'n' (number of rows to display) was wrong!", call. = FALSE )
				}
				# n is numeric and in correct range
			} else if( n != "all" ){
				stop( "Could not understand the choice of 'n'", call. = FALSE )
			} else {
				# n == "all"
				n <- nrow( cov.data )
			}
		}

		if( !missing( from ) ){
			from.given <- TRUE
			if( !is.numeric( from ) ){
				stop( "'From' should be numeric!", call. = FALSE )
			}
		}

		if( !missing( to ) ){
			if( !is.numeric( to ) ){
				stop( "'To' should be numeric!", call. = FALSE )
			}
			if( n.given & from.given ){
				warning( "You specified too many parameters! Chosing only 'from' and 'to'." )
			} else if( n.given ) {
				# only 'n' and 'to' are specified
				from <- to - n + 1
			} else if( !from.given ) {
				# only 'to' is specified
				from <- 1
			}
		} else if( n.given & from.given ) {
		# only 'from' and 'n' is specified
			to <- from + n - 1
		} else if( from.given ) {
		# only 'from' is specified
			to <- nrow( cov.data )
		} else if( n.given ) {
		# only 'n' is specified
			from <- 1
			to <- from + n - 1
		} else {
			# nothing specified: print the first 5 entries
			from <- 1
			to <- n
		}

		if( from < 1 | to > nrow( cov.data ) | from > to ){
			stop( paste0( "Wrong specification of 'from' (", from ,") and 'to' (", to, ")!" ), call. = FALSE )
		}
		data.out <- cov.data[ from:to, ]
	}
	return( data.out )
}
