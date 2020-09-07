f.insert.row <- function( df, newrow, r ){
	if( r != ( nrow( df ) + 1 ) ){
		df[ seq( r + 1, nrow( df ) + 1 ), ] <- df[ seq( r, nrow( df ) ), ]
	}
	df[ r, ] <- newrow
	return( df )
}
