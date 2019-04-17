f.make.out.filename <- function( file.in, file.out, dir.out, root = "gen", overwrite = NULL ){
	if( is.null( file.out ) ){
		fname.vec <- unlist( strsplit( file.in, split = "/", fixed = TRUE ) )
		l.fname.vec <- length( fname.vec )
		if( l.fname.vec > 1 ){
			file.in <- fname.vec[ l.fname.vec ]
		}
		fname.vec <- unlist( strsplit( file.in, split = ".", fixed = TRUE ) )
		l.fname.vec <- length( fname.vec )
		if( l.fname.vec > 1 ){
			file.out <- paste( fname.vec[ -l.fname.vec ], collapse = "_" )
		} else {
			file.out <- file.in
		}
	}
	
	file.out.base <- paste0( file.out, "_", root )
	file.out.ff <- paste0( dir.out, "/", file.out, "_", root, ".ffData")
	file.out.aux <- paste0( dir.out, "/", file.out, "_", root, ".RData" )
	if( file.exists( file.out.ff ) | file.exists( file.out.aux ) ){
		cat( "The output file(s) exist! \n" )
		if( is.null( overwrite ) ){
			answer <- readline( paste( "Do you want to overwrite file(s)? (y/n)" ) )
			if( tolower( answer ) != "y" ){
				stop( "Stopped without overwriting files.", call. = FALSE )
				return( NULL )
			}
		} else if( !is.logical( overwrite ) | !overwrite ){
			stop( "Stopped without overwriting files.", call. = FALSE )
			return( NULL )
		} else if( overwrite ){
			ff::ffdrop( paste0( dir.out, "/", file.out.base ) )
		}
	}
	
	return( list( file.out.base = file.out.base, file.out.ff = file.out.ff, file.out.aux = file.out.aux ) )
}
