#' Getter for the list with genetic or environmental data
#'
#' This function returns the columns requested from the column-wise splitted genetic/environmental data in ff format.
#'
#' @param gen.data The list with genetic/environmental data (in ff format).
#' @param cols The sequence numbering (or listing) the columns which are to be retrieved; this should enumerate the columns as if all the data was in one large matrix.
#' @param by.colname Logical: is 'cols' given as a number sequence (FALSE, default) or character vector (TRUE) naming the columns to choose?
#'
#' @return Matrix (in ff format) with all the requested columns.
#'
#' @keywords internal

f.get.gen.data.cols <- function( gen.data, cols, by.colname = FALSE ){
	ncols.per.chunk <- ncol( gen.data[[ 1 ]] )
	if( by.colname ){
		gen.list.info <- f.get.which.gen.el.names( cols, gen.data )
	} else {
		gen.list.info <- f.get.which.gen.el( cols, ncols.per.chunk )
	}
	which.gen.chunk <- gen.list.info$chunk.no
	which.cols.chunk <- gen.list.info$col.no
	
	new.dim <- c( nrow( gen.data[[ 1 ]] ), length( cols ) )

	all.levels <- c()
	for( i in unique( which.gen.chunk ) ){
		all.levels <- union( all.levels, ff::levels.ff( gen.data[[ i ]] ) )
	}
	
	new.colnames <- c()
	gen.data.out <- ff::ff( NA, levels = all.levels, dim = new.dim, vmode = ff::vmode( gen.data[[ 1 ]] ) )
	for( i in 1:length( cols ) ){
		gen.data.out[ ,i ] <- gen.data[[ which.gen.chunk[ i ] ]][ ,which.cols.chunk[ i ] ]
		new.colnames <- c( new.colnames, 
			colnames( gen.data[[ which.gen.chunk[ i ] ]] )[ which.cols.chunk[ i ] ] )
	}
	colnames( gen.data.out ) <- new.colnames
	return( gen.data.out )
}
