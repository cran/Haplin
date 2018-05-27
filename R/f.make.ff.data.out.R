make.ff.data.out <- function( covd, gend, covd.names, gend.names, data.as.is = FALSE ){
	if( !data.as.is ){
		if( !missing( covd ) ){
			covd <- dframe( covd )
		}
		
		gend.char.mat <- as.matrix( gend, mode = "character" )
		gen.levels <- unique( as.character( gend.char.mat ) )
		if( !any( is.na( gen.levels ) ) ){
			gen.levels <- c( gen.levels, NA )
		}
		gend <- ff::as.ff( gend.char.mat, vmode = .haplinEnv$.vmode.gen.data, levels = gen.levels )
	}

	if( !missing( covd ) ){
		data <- list( cov.data = covd, gen.data = gend, aux = list() )
	} else {
		data <- list( cov.data = NULL, gen.data = gend, aux = list() )
	}

	if( !missing( covd.names ) ){
		colnames( data$cov.data ) <- covd.names
	}
	if( !missing( gend.names ) ){
		colnames( data$gen.data ) <- gend.names
	}

#class(data) <- "haplin.ready" # xxx tja
	return( data )
}
