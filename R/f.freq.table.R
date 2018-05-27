f.freq.table <- function( x, withNA = FALSE ){
	.tab <- table( as.character( x ), useNA = "always" )
	.nas <- .tab[ is.na( names( .tab ) ) ]
	.tab <- .tab[ .tab != 0 ]
	if( !withNA ){
		.tab <- .tab[ !is.na( names( .tab ) ) ]
	}
	# sorting by alphabet
	sort.names <- sort( names( .tab ), index.return = TRUE )
	no.na.ids <- which( !is.na( names( .tab ) ) )
	if( withNA ){
		.tab <- .tab[ c( no.na.ids[ sort.names$ix ], which( is.na( names( .tab ) ) ) ) ]
	} else {
		.tab <- .tab[ no.na.ids[ sort.names$ix ] ]
	}
	return( list( tab = .tab, nas = .nas ) )
}
