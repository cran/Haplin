#' Extract a part of an ff array (matrix) and return it as numeric matrix
#'
#' This function returns the requested columns and/or rows and converts it to a numeric matrix (by default, the extraction from an ff matrix with defined levels returns a factor vector).
#'
#' @param gen.data The ff matrix with genetic data.
#' @param cols The sequence numbering the columns which are to be retrieved, if NULL (default) then all the columns will be extracted.
#' @param rows The sequence numbering the rows which are to be retrieved, if NULL (default) then all the rows will be extracted.
#'
#' @return Matrix with all the requested columns, data in numeric representation.
#'
#' @keywords internal

f.extract.ff.numeric <- function( gen.data, cols = NULL, rows = NULL ){
	if( is.null( dim( gen.data ) ) ){
		stop( "'gen.data' should be a matrix!" )
	}

	if( is.null( cols ) & is.null( rows ) ){
		.data.extr.num <- as.numeric( levels( gen.data[,] ) )[ gen.data[,] ]
		out.cols <- ncol( gen.data )
	} else {
		if( !is.null( cols ) ){
			if( !all( cols %in% 1:ncol( gen.data ) ) ){
				stop( "Some requested columns are outside the column range of gen.data!" )
			}
			
			if( !is.null( rows ) ){
				if( !all( rows %in% 1:nrow( gen.data ) ) ){
				stop( "Some requested rows are outside the row range of gen.data!" )
				}
				.data.extr <- gen.data[ rows,cols ]
			} else {
				.data.extr <- gen.data[ ,cols ]
			}
			out.cols <- length( cols )
		} else { # cols = NULL
			if( !all( rows %in% 1:nrow( gen.data ) ) ){
				stop( "Some requested rows are outside the row range of gen.data!" )
			}
			.data.extr <- gen.data[ rows, ]
			out.cols <- ncol( gen.data )
		}
		.data.extr.num <- as.numeric( levels( .data.extr ) )[ .data.extr ]
	}
	
	return( matrix( .data.extr.num, ncol = out.cols ) )
}
