f.ped.to.mfc.new <- function( data.in, design ){
	gen.list.length <- length( data.in$gen.data )

	# create unique ids for individuals
	cat( "Converting PED format to internal haplin...\n" )
	
	# if the children codes are not unique
	all.id.c <- table( data.in$cov.data[ ,"id.c" ] )
	if( any( all.id.c > 1 ) ){
		cat( "   Creating unique IDs for individuals...\n" )
		orig.cov.colnames <- colnames( data.in$cov.data )
		data.in$cov.data <- t( apply( data.in$cov.data, 1, function( x ){
			if( x[ 4 ] == 0 | x[ 3 ] == 0 ){
				new.ids <- c( paste( x[ 1 ], x[ 2 ], sep = "_" ), x[ 3:4 ] )
			} else {
				new.ids <- paste( x[ 1 ], x[ 2:4 ], sep = "_" )
			}
			return( c( x[ 1 ], new.ids, x[ 5:length( x ) ] ) )
		} ) )
		colnames( data.in$cov.data ) <- orig.cov.colnames
		cat( "   ...done.\n" )
	}
	
	id <- data.in$cov.data[ ,"id.c" ]
	# sort the families and check coding
	cat( "   Sorting and re-coding families...\n" )
	pedIndex <- f.prep.pedIndex( data.in$cov.data )

	if( design %in% c( "triad", "cc.triad" ) ){
		# SELECT LINES OF data CORRESPONDING TO EITHER MOTHER, FATHER OR CHILD
		# NOTE THAT DATA LINES NOT CORRESPONDING TO INDIVIDUALS IDENTIFIED IN THE
		# pedIndex FILE WILL NOT BE SELECTED.
		d.m <- match( pedIndex[ ,'id.mother' ], id )
		d.f <- match( pedIndex[ ,'id.father' ], id )
		d.c <- match( pedIndex[ ,'id.child' ], id )
		
		gen.data.ut <- list()
		for( i in 1:gen.list.length ){
			cur.gen.cols <- data.in$gen.data[[ i ]]
			# DIMENSION FOR NEW DATA SET
			new.dim <- c( nrow( pedIndex ), 3 * ncol( cur.gen.cols ) )
			# NEW INTERLACED DATA SET
			all.levels <- levels( cur.gen.cols[,] )
			
			cur.gen.data.ut <- ff::ff( dim = new.dim, vmode = .haplinEnv$.vmode.gen.data, levels = all.levels )
			seq.m <- c( rbind( seq(1, new.dim[2], 6), seq(2, new.dim[2], 6) ) )
			if( any( is.na( d.m ) ) ){
				which.m.na <- which( is.na( d.m ) )
				tmp.genes <- data.frame( matrix( as.character( cur.gen.cols[ d.m[ -which.m.na ], ] ), nrow = length( d.m[ -which.m.na ] ) ) )
				na.row <- rep( NA, ncol( tmp.genes ) )
				for( l in which.m.na ){
					tmp.genes <- f.insert.row( tmp.genes, na.row, l )
				}
				tmp.genes <- unlist( tmp.genes, use.names = FALSE )
			} else {
				tmp.genes <- cur.gen.cols[ d.m, ]
			}
			cur.gen.data.ut[ ,seq.m ] <- tmp.genes
			seq.f <- c( rbind( seq(3, new.dim[2], 6), seq(4, new.dim[2], 6) ) )
			if( any( is.na( d.f ) ) ){
				which.f.na <- which( is.na( d.f ) )
				tmp.genes <- data.frame( matrix( as.character( cur.gen.cols[ d.f[ -which.f.na ], ] ), nrow = length( d.f[ -which.f.na ] ) ) )
				na.row <- rep( NA, ncol( tmp.genes ) )
				for( l in which.f.na ){
					tmp.genes <- f.insert.row( tmp.genes, na.row, l )
				}
				tmp.genes <- unlist( tmp.genes, use.names = FALSE )
			} else {
				tmp.genes <- cur.gen.cols[ d.f, ]
			}
				cur.gen.data.ut[,seq.f] <- tmp.genes
			seq.c <- c( rbind( seq(5, new.dim[2], 6), seq(6, new.dim[2], 6) ) )
			cur.gen.data.ut[,seq.c] <- cur.gen.cols[ d.c, ]
			gen.data.ut <- c( gen.data.ut, list( cur.gen.data.ut ) )
		}
		
		# now - the covariate data
		new.dim <- c( nrow( pedIndex ), 3 * ( ncol( data.in$cov.data ) - 4 ) + 4 )
		cov.data.reordered <- matrix( NA_character_, nrow = new.dim[1], ncol = new.dim[2] )
		chosen.cols <- c( "id.c",colnames( data.in$cov.data )[-(1:4)] )
		cov.data.reordered[,1] <- pedIndex[ ,"family" ]
		cov.data.reordered[, seq(2, new.dim[2], 3)] <- as.matrix( data.in$cov.data[ d.m,chosen.cols ], mode = "character" )
		cov.data.reordered[, seq(3, new.dim[2], 3)] <- as.matrix( data.in$cov.data[ d.f,chosen.cols ], mode = "character" )
		cov.data.reordered[, seq(4, new.dim[2], 3)] <- as.matrix( data.in$cov.data[ d.c,chosen.cols ], mode = "character" )
		
		labs <- c("m", "f", "c")
	} else if( design == "cc" ){
		# SELECT LINES OF data CORRESPONDING TO CHILDREN IN pedIndex FILE
		# NOTE THAT DATA LINES NOT CORRESPONDING TO INDIVIDUALS IDENTIFIED IN THE
		# pedIndex FILE WILL NOT BE SELECTED.
		d.c <- match( pedIndex[ ,'id.child' ], id )
		gen.data.ut <- list()
		for( i in 1:gen.list.length ){
			cur.gen.cols <- data.in$gen.data[[ i ]]
			# DIMENSION FOR NEW DATA SET
			new.dim <- c( nrow( pedIndex ), ncol( cur.gen.cols ) )
			# NEW INTERLACED DATA SET
			all.levels <- levels( cur.gen.cols[,] )

			cur.gen.data.ut <- ff::ff( cur.gen.cols[ d.c, ], dim = new.dim, vmode = .haplinEnv$.vmode.gen.data, levels = all.levels )
			gen.data.ut <- c( gen.data.ut, list( cur.gen.data.ut ) )
		}
		
		chosen.cols <- c( "id.c",colnames( data.in$cov.data )[-(1:4)] )
		cov.data.reordered <- as.matrix( data.in$cov.data[ d.c,chosen.cols ], mode = "character" )
		cov.data.reordered <- cbind( pedIndex[ ,"family" ], cov.data.reordered )
		
		labs <- "c"
	}
	
	for( i in 1:gen.list.length ){
		cur.gen.cols <- data.in$gen.data[[ i ]]
		
		markers1 <- grep( "_a", colnames( cur.gen.cols ) )
		markers2 <- grep( "_b", colnames( cur.gen.cols ) )

		marker.names.a <- as.vector( t( outer( colnames( cur.gen.cols )[ markers1 ], labs, paste, sep = "_" ) ) )
		marker.names.b <- as.vector( t( outer( colnames( cur.gen.cols )[ markers2 ], labs, paste, sep = "_" ) ) )
		gen.data.colnames <- as.vector( rbind( marker.names.a, marker.names.b ) )

		colnames( gen.data.ut[[ i ]] ) <- gen.data.colnames
	}

	n.vars <- ncol( data.in$cov.data )
	orig.cov.colnames <- colnames( data.in$cov.data )
	cov.data.colnames <- c( paste( "id", labs, sep = "." ), paste( "sex", labs, sep = "." ), paste( "cc", labs, sep = "." ) )
	if( n.vars > 6 ){
		cov.data.colnames <- c( cov.data.colnames, sapply( orig.cov.colnames[ -(1:6) ], function( x ){ 
			paste( x, labs, sep = "." )
		} ) )
	}
	colnames( cov.data.reordered ) <- c( "id.fam", cov.data.colnames )

	return( list( cov.data = cov.data.reordered, gen.data = gen.data.ut ) )
}
