#' Display chosen genotypes
#'
#' This is a help function to extract genotypes from an object read in with
#'  \link{genDataRead} (or loaded with \link{genDataLoad}).
#'
#' @param data.in The data read in by \link{genDataRead}.
#' @param design The design used in the study - choose from:
#'   \itemize{
#'     \item \emph{triad} - data includes genotypes of mother, father and child;
#'     \item \emph{cc} - classical case-control;
#'     \item \emph{cc.triad} - hybrid design: triads with cases and controls;
#'   }.
#' @param n Number of rows to display or "all" (default: 5).
#' @param from From which row to display (optional, default: from the first).
#' @param to To which row to display (optional).
#' @param sex If the sex column is part of the phenotypic information, the user can
#'    choose based on one of the categories used in this column (optional);
#'    NB: this does not combine with the 'to' and 'from' arguments.
#' @param markers A vector specifying which markers to display or "all" (default: first
#'    5); NB: the user can specify the markers by numbers or by their names.
#'
#' @return A table with genotypes extracted from 'data.in'.
#'

### @param ids Individual IDs to display (optional); NB: this does not combine with the
###    'to' and 'from' arguments.

showGen <- function( data.in, design = "triad", n = 5, from, to, sex, markers = 1:5 ){
	# check if input data is in correct format
	if( !is( data.in, "haplin.data" ) ||
	  !all( names( data.in ) == .haplinEnv$.haplin.data.names ) ){
		stop( "The input data is not in the correct format!", call. = FALSE )
	}
	
	cov.data <- data.in$cov.data
	gen.data <- data.in$gen.data
	n.given <- FALSE
	from.given <- FALSE
	which.rows.show <- TRUE

	if( !missing( sex ) ){
		# check if there is any phenotypic info
		if( is.null( cov.data ) ){
			stop( "There is no phenotypic/covariate information in your data!", call. = FALSE )
		}
		if( !( "sex" %in% colnames( cov.data ) ) ){
			stop( "There is no 'sex' column in the dataset!" )
		}
		if( !( sex %in% as.numeric( unique( cov.data[ ,"sex" ] ) ) ) ){
			stop( paste( "Wrong sex chosen - 'sex' column has the following categories:", paste( sort( unique( cov.data[ ,"sex" ] ) ), collapse = "," ) ), call. = FALSE )
		}
		
		which.rows.show <- as.numeric( cov.data[ ,"sex" ] ) == sex
	} else {
		if( !missing( n ) ){
			n.given <- TRUE
			if( is.numeric( n ) ){
				if( n < 1 | n > nrow( gen.data[[ 1 ]] ) ){
					stop( "Your specification of 'n' (number of rows to display) was wrong!", call. = FALSE )
				}
				# n is numeric and in correct range
			} else if( n != "all" ){
				stop( "Could not understand the choice of 'n'", call. = FALSE )
			} else {
				# n == "all"
				n <- nrow( gen.data[[ 1 ]] )
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
			to <- nrow( gen.data[[ 1 ]] )
		} else if( n.given ) {
		# only 'n' is specified
			from <- 1
			to <- from + n - 1
		} else {
			# nothing specified: print the first 5 entries
			from <- 1
			to <- min( n, nrow( gen.data[[ 1 ]] ) )
		}

		if( from < 1 | to > nrow( gen.data[[ 1 ]] ) | from > to ){
			stop( paste0( "Wrong specification of 'from' (", from ,") and 'to' (", to, ")!" ), call. = FALSE )
		}
		which.rows.show <- from:to
	}
	
	format <- data.in$aux$info$filespecs$format
	family <- "c"
	ncols.per.marker <- 2
	if( design %in% c( "triad", "cc.triad" ) & format == "haplin" ){
		family <- "mfc"
		ncols.per.marker <- 6
	}

	max.markers <- nsnps( data.in )
	all.markers <- max.markers*ncols.per.marker
	if( is.numeric( markers ) ){
		markers.by.name <- FALSE
		if( max( markers ) > max.markers ){
			warning( "It appears that your data has less markers (", max.markers,") than requested (", max( markers ), "), adjusting accordingly." )
			orig.markers <- markers
			which.markers.higher.max <- orig.markers > max.markers
			markers <- orig.markers[ !which.markers.higher.max ]
		}

		all.cols.names <- unlist(lapply(gen.data, colnames))
		
		markers <- f.sel.markers(
		  n.vars = 0,
		  markers = markers,
		  family = family,
		  split = TRUE,
		  all.marker.names = all.cols.names
		)
	} else if( length( markers ) == 1 & markers == "all" ){
		markers.by.name <- FALSE
		if( length( gen.data ) > 1 ){
			answer <- readline( paste( "You are requesting a lot of data. Are you sure? (y/n)" ) )
			if( tolower( answer ) != "y" ){
				stop( "Stopped as per request.", call. = FALSE )
			}
		}
		
		markers <- 1:all.markers
	} else {
		markers.by.name <- TRUE
		stop( "Not implemented yet!" )
		# TODO: change the f.sel.markers function?
	}
	data.out <- f.get.gen.data.cols( gen.data, cols = markers, by.colname = markers.by.name )
	data.out <- as.ff( data.out[ which.rows.show, ], vmode = .haplinEnv$.vmode.gen.data )
	
	return( data.out )
}
