#' Checking uniqueness of individuals' IDs
#'
#' This checks whether the individuals have unique IDs. If not, it creates new IDs by
#'   adding the family ID to the indiv. ID. Since the function checks for specific
#'   column names, this basically means that the input should be as read from 
#'   a PED file.
#'
#' @param data.cov A matrix with "id.fam", "id.c", "id.m", and "id.f" columns.
#'
#' @return List with:
#'   \itemize{
#'     \item ids - new IDs in the same order as "id.c" column in the input data
#'     \item pedIndex - matrix with "family", "id.child", "id.mother", and
#'        "id.father" columns
#'   }
#'
#' @keywords internal

f.check.unique.ids <- function( data.cov ){
	all.id.c <- table( data.cov[ ,"id.c" ] )
	if( any( all.id.c > 1 ) ){
		cat( "   Creating unique IDs for individuals...\n" )
		orig.cov.colnames <- colnames( data.cov )
		data.cov <- t( apply( data.cov, 1, function( x ){
			if( x[ 4 ] == 0 | x[ 3 ] == 0 ){
				new.ids <- c( paste( x[ 1 ], x[ 2 ], sep = "_" ), x[ 3:4 ] )
			} else {
				new.ids <- paste( x[ 1 ], x[ 2:4 ], sep = "_" )
			}
			return( c( x[ 1 ], new.ids, x[ 5:length( x ) ] ) )
		} ) )
		colnames( data.cov ) <- orig.cov.colnames
		cat( "   ...done.\n" )
	}
	
	id <- data.cov[ ,"id.c" ]
	# sort the families and check coding
	pedIndex <- f.prep.pedIndex( data.cov )

	return( list( ids = id, pedIndex = pedIndex ) )
}
