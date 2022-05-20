#' Intern function for creating column names for genotype data
#' 
#' Creating column names either from a .map file or generating dummy names
#'
#' If the .map file is given, the SNP names are read in, if not dummy names are
#'  created by attaching numbers to "m" prefix. Then, a "_a" and "_b" suffix is
#'  attached to each allele, respectively. Finally, if the design is "triad" or
#'  "cc.triad", the following suffixes are attached: "_m" for the mother's
#'  alleles, "_f" for the father's, and "_c" for the child's.
#'
#' @param map.file Filename (with path if the file is not in current directory) of the
#'   .map file holding the SNP names, if available (see Details).
#' @param ncol Number of columns IN TOTAL in the dataset containing only the genotype data
#' @param format Format of data (will influence how data is processed) - choose from:
#'   \itemize{
#'     \item \emph{haplin} - data already in one row per family,
#'     \item \emph{ped} - data from .ped file, each row represents an individual.
#'   }.
#' @param design The design used in the study - choose from:
#'   \itemize{
#'     \item \emph{triad} - (default), data includes genotypes of mother, father and child;
#'     \item \emph{cc} - classical case-control;
#'     \item \emph{cc.triad} - hybrid design: triads with cases and controls
#'   }.
#'
#' @section Details:
#' The .map file should contain at least two columns, where the second one contains SNP 
#'   names. Any additional columns should be separated by a whitespace character, but will 
#'   be ignored. The file should contain a header.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \emph{gen.data.colnames} - a vector with names of columns, length equal to 
#'       the number of columns in the genotype dataset (i.e., depending on the format and
#'       design).
#'     \item \emph{marker.names} - a vector containing the names of markers, as read in
#'       from 'map.file', or dummy names.
#'   }
#'

f.create.snp.names <- function( map.file, ncol, format, design ){
	cat( "Reading the marker names... \n" )

	if( format == "ped" ){
		ncol.per.locus <- 2
	} else if( design %in% c( "triad", "cc.triad" ) ){
		ncol.per.locus <- 6
	} else {
		ncol.per.locus <- 2
	}

	marker.names <- c()
	if( !is.null( map.file ) ){
		marker.names <- read.table( map.file, header = TRUE, stringsAsFactors = FALSE )
		if( ( nrow( marker.names ) * ncol.per.locus ) != ncol ){
			marker.names <- c()
		}
	}
	if( length( marker.names ) == 0 ){
		warning( "No map file given, map file empty or the number of map file rows not equal to the number of markers in data; will generate dummy marker names.", call. = FALSE )
		marker.names <- paste( "m", 1:( ncol/ncol.per.locus ), sep = "" )
	} else {
		marker.names <- marker.names[ ,2 ]
	}

	# all the marker names will start with l_
	gen.data.colnames <- paste( "l", as.character( marker.names ), sep = "_" )
	markers1 <- paste( gen.data.colnames, "a", sep = "_" )
	markers2 <- paste( gen.data.colnames, "b", sep = "_" )

	if( format == "ped" | ( format == "haplin" & design == "cc" ) ){
		gen.data.colnames <- as.vector( rbind( markers1, markers2 ) )
	} else if( format == "haplin" & design %in% c( "triad", "cc.triad" ) ){
		labs <- c("m", "f", "c")
		marker.names.a <- as.vector( t( outer( markers1, labs, paste, sep = "_" ) ) )
		marker.names.b <- as.vector( t( outer( markers2, labs, paste, sep = "_" ) ) )
		gen.data.colnames <- as.vector( rbind( marker.names.a, marker.names.b ) )
	} else {
		stop( "Problem with design and format!", call. = FALSE )
	}
	cat( "...done\n" )
	
	return( list( gen.data.colnames = gen.data.colnames, marker.names = marker.names ) )
}
