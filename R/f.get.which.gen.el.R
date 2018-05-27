f.get.which.gen.el <- function( cols, ncols.per.chunk ){
	if( !is.numeric( cols ) ){
		stop( "The 'cols' argument should be numeric!" )
	}

	which.gen.chunk <- vector( length( cols ), mode = "integer" )
	which.cols.chunk <- vector( length( cols ), mode = "integer" )
	for( i in 1:length( cols ) ){
		which.gen.chunk[ i ] <- ( cols[ i ] - 1 ) %/% ncols.per.chunk + 1
		which.cols.chunk[ i ] <- ( cols[ i ] - ( which.gen.chunk[ i ] - 1 )*ncols.per.chunk )%%( ncols.per.chunk + 1 )
	}
	
	return( list( chunk.no = which.gen.chunk, col.no = which.cols.chunk ) )
}
