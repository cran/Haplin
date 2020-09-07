#' Count the number of individuals in the data
#'
#' This is a help function to count the number of individuals in an object read in with
#'  \link{genDataRead} (or loaded with \link{genDataLoad}).
#'
#' @param data.in The data read in by \link{genDataRead}.
#' @param design The design used in the study - choose from:
#'   \itemize{
#'     \item \emph{triad} (default) - data includes genotypes of mother, father and child;
#'     \item \emph{cc} - classical case-control;
#'     \item \emph{cc.triad} - hybrid design: triads with cases and controls
#'   }
#'
#' @return How many individuals (integer).
#'
nindiv <- function( data.in, design = "triad" ){
	# check if input data is in correct format
	if( !is( data.in, "haplin.data" ) ||
	  !all( names( data.in ) == .haplinEnv$.haplin.data.names ) ){
		stop( "The input data is not in the correct format!", call. = FALSE )
	}

	design.list <- get( ".design.list", envir = .haplinEnv )
	if( !( design %in% design.list ) ){
		stop( "Given design(", design,") not recognized! Design has to be one of: ", paste( design.list, collapse = ", " ), call. = FALSE )
	}

	tot.no.lines <- nrow( data.in$gen.data[[1]] )
	format <- data.in$aux$info$filespecs$format
	nindiv.per.line <- 1 # true only for PED format and haplin cc
	if( design %in% c( "triad", "cc.triad" ) & format == "haplin" ){
		nindiv.per.line <- 3
	}
	return( tot.no.lines * nindiv.per.line )
}
