#' Converting Haplin-formatted matrix into the new format
#'
#' Internal function for converting an R matrix object into the new format of data used in Haplin. The output is the same as from \link{genDataRead}.
#'
#' @param data A character matrix (NB: it's assumed that this matrix contains genotypes
#'   and covariate data, if any; i.e., as if it was a table read in from a haplin-formatted
#'   file).
#' @param n.vars Number of columns with covariate variables - if the data does not contain
#'   covariates, give 0 explicitly.
#' @param cov.header Optional; if there are covariates in the data, you can give their 
#'   names explicitly here. Otherwise, dummy names will be created.
#' @param gen.levels Optional; a vector with all the possible values for alleles in the
#'   genotype part of data. If not given, these will be assessed from the given data.
#'
#' @return A list object with two elements:
#'   \itemize{
#'     \item \emph{cov.data} - a \code{data.frame} with covariate data (if available in
#'        the input file)
#'     \item \emph{gen.data} - a list with chunks of the genetic data; the data is divided
#'        column-wise, using 10,000 columns per chunk; each element of this list is a
#'        \link[ff]{ff} matrix
#'   }
#'
f.convert.matrix.ff <- function( data = stop( "You must give the data to convert!", call. = FALSE ), n.vars = stop( "You must explicitly give the number of columns with covariates!", call. = FALSE ), cov.header, gen.levels ){
	if( n.vars == 0 ){
		gen.data <- data
	} else {
		cov.data <- data[ ,1:n.vars ]
		if( !missing( cov.header ) ){
			if( length( cov.header ) != ncol( cov.data ) ){
				stop( "You gave the header for covariate data but its length doesn't match the number of covariate columns!", call. = FALSE )
			}
			colnames( cov.data ) <- cov.header
		} else {
			cat( "Will generate dummy names for covariates.\n" )
			colnames( cov.data ) <- paste0( "cov.", 1:n.vars )
		}
		gen.data <- data[ ,-(1:n.vars) ]
	}
	
	if( missing( gen.levels ) ){
		gen.levels <- unique( as.character( gen.data ) )
	}
	if( !any( is.na( gen.levels ) ) ){
		gen.levels <- c( gen.levels, NA )
	}
	
	nb.cols.per.chunk <- get( ".nb.cols.per.chunk", envir = .haplinEnv )
	nb.cols.gen.data <- ncol( gen.data )
	nb.col.chunks <- ceiling( nb.cols.gen.data / nb.cols.per.chunk )
	nb.rows.tot <- nrow( gen.data )
	
	gen.data.col.wise <- list()
	for( i in 1:nb.col.chunks ){
		cur.cols <- ( ( i-1 )*nb.cols.per.chunk + 1 ):( min( i*nb.cols.per.chunk, nb.cols.gen.data ) )
		tmp.gen.data <- ff::ff( gen.data[ ,cur.cols ], vmode = .haplinEnv$.vmode.gen.data, levels = gen.levels, dim = c( nb.rows.tot, min( nb.cols.per.chunk, max( cur.cols ) - min( cur.cols ) + 1 ) ) )
		gen.data.col.wise <- c( gen.data.col.wise, list( tmp.gen.data ) )
		rm( tmp.gen.data )
	}
	
	if( n.vars == 0 ){
		out.data <- make.ff.data.out( gend = gen.data.col.wise, data.as.is = TRUE )
	} else {
		out.data <- make.ff.data.out( covd = cov.data, gend = gen.data.col.wise, data.as.is = TRUE )
	}
	class( out.data ) <- "haplin.data"
	
	return( out.data )
}
