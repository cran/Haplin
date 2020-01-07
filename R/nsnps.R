#' Count the number of markers in the data
#'
#' This is a help function to count the number of markers in an object read in with
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
#' @return How many markers (integer).
#'
nsnps <- function( data.in, design = "triad" ){
	# check if input data is in correct format
	if( !is( data.in, "haplin.data" ) ||
	  !all( names( data.in ) == .haplinEnv$.haplin.data.names ) ){
		stop( "The input data is not in the correct format!", call. = FALSE )
	}
	
	design.list <- get( ".design.list", envir = .haplinEnv )
	if( !( design %in% design.list ) ){
		stop( "Given design(", design,") not recognized! Design has to be one of: ", paste( design.list, collapse = ", " ), call. = FALSE )
	}
	
	ncols.per.chunk <- sapply( data.in$gen.data, ncol )
	
	format <- data.in$aux$info$filespecs$format
	ncols.per.marker <- 2
	if( design %in% c( "triad", "cc.triad" ) & format == "haplin" ){
		ncols.per.marker <- 6
	}
	return( sum( ncols.per.chunk ) / ncols.per.marker )
}
