#' Getter for all rows with fathers' data
#'
#' Wrapper function for \link{genDataGetPart} that returns a subset of the data containing 
#'  only fathers.
#'
#' @param data.in The data object (in format as the output of \link{genDataRead}); note
#'   that the design of the data is assumed to be triad.
#' @param file.out The base for the output filename (default: "my_data_onlyFathers").
#' @param dir.out The path to the directory where the output files will be saved.
#' @param overwrite Whether to overwrite the output files: if NULL (default), will prompt
#'   the user to give answer; set to TRUE, will automatically overwrite any existing files;
#'   and set to FALSE, will stop if the output files exist.
#'
#' @return A list object with three elements:
#'   \itemize{
#'     \item \emph{cov.data} - a \code{data.frame} with covariate data (if available in
#'        the input file)
#'     \item \emph{gen.data} - a list with chunks of the genetic data; the data is divided
#'        column-wise, using 10,000 columns per chunk; each element of this list is a
#'        \link[ff]{ff} matrix
#'     \item \emph{aux} - a list with meta-data and important parameters.
#'   }
#'   This now contains only the selected subset of data.
#'

getFathers <- function( data.in = stop( "No data given!", call. = FALSE ), file.out = "my_data_onlyFathers", dir.out = ".", overwrite = NULL ){
	format <- data.in$aux$info$filespecs$format
	
	if( format == "ped" ){
		# first - get all the IDs and check family structure
		new.ids <- f.check.unique.ids( data.in$cov.data )
		id <- new.ids$ids
		pedIndex <- new.ids$pedIndex
		
		# which IDs are fathers?
		d.c <- match( pedIndex[ ,'id.father' ], id )
		return( genDataGetPart( data.in, design = "triad", rows = d.c, file.out = file.out, dir.out = dir.out, overwrite = overwrite ) )
	} else if( format == "haplin" ){
		stop( "Not implemented yet", call. = FALSE )
	} else {
		stop( paste( "Unrecognized format:", format ), call. = FALSE )
	}
}
