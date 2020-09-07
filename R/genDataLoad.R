#' Loading the data previously read in and saved by "genDataRead"
#'
#' This function loads the data from the saved .ffData and .RData files,
#' and prepares the data to subsequent analysis.
#'
#' @param filename The base of the filenames; i.e. if the data is saved in
#'  "my_data_gen.ffData", "my_data_gen.RData" and "my_data_cov.RData", then the 'filename'
#'  should be "my_data".
#' @param dir.in The path to the directory where files were saved (defaults to the
#'  current directory).
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
#'
genDataLoad <- function( filename = stop( "'filename' must be given!" ), dir.in = "." ){
	file.in.ff <- file.path( dir.in, paste0( filename, "_gen.ffData" ) )
	if( !file.exists( file.in.ff ) ){
		stop( "The file(s) doesn't seem to exist!", call. = FALSE )
	}
	
	gen.cols.name <- get( ".gen.cols.name", envir = .haplinEnv )
	file.in.base <- file.path( dir.in, paste0( filename, "_gen" ) )
	
	suppressWarnings( ff::ffload( file.in.base, rootpath = getOption( "fftempdir" ) ) )
	loaded.objects <- ls( pattern = paste0( gen.cols.name, ".[[:digit:]]" ) )
	loaded.aux <- ls( pattern = "aux" )
	if( length( loaded.aux ) == 0 ){
		stop( "Problem: no 'aux'!" )
	}
	aux <- get( "aux" )
	
	# maintain the correct order!
	if( length( loaded.objects ) > 9 ){
		chunks.no <- sapply( loaded.objects, function(x){
			tmp.string <- unlist( strsplit( x, split = ".", fixed = TRUE ) )
			return( as.numeric( tmp.string[ length( tmp.string ) ] ) )
		}, USE.NAMES = FALSE )
		sorted.chunks <- sort( chunks.no, index.return = TRUE )
		loaded.objects <- loaded.objects[ sorted.chunks$ix ]
	}
	
	gen.data.col.wise <- list()
	for( i in loaded.objects ){
		cur.chunk <- get( i )
		if( !is.null( dim( cur.chunk ) ) ){
			gen.data.col.wise <- c( gen.data.col.wise, list( cur.chunk ) )
		}
		rm( cur.chunk )
	}
	
	if( !exists( "cov.data.in" ) ){
		cov.data.in <- NULL
	}

	data.out <- list( cov.data = cov.data.in, gen.data = gen.data.col.wise, aux = aux )
	class( data.out ) <- aux$class
	
	return( data.out )
}
