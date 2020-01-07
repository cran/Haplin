f.get.which.gen.el.names <- function( cols, gen.data ){
	if( !is.character( cols ) ){
		stop( "The 'cols' argument should be a character vector!" )
	}

	out.list <- sapply( cols, function( col.name ){
		for( chunk.no in 1:length( gen.data )){
			which.col <- match( col.name, colnames( gen.data[ chunk.no ] ), nomatch = -1 )
			if( which.col != -1 ){
				return( list( chunk.no = chunk.no, col.no = which.col ) )
			}
		}
	} )
	
	out.df <- data.frame( chunk.no = as.numeric( out.list[ "chunk.no", ] ), col.no = as.numeric( out.list[ "col.no", ] ) )
	return( out.df )
}
