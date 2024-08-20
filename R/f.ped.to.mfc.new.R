f.ped.to.mfc.new <- function( data.in, design ){
	cat( "Converting PED format to internal haplin...\n" )
	#
	gen.list.length <- length( data.in$gen.data )
	#
	# if the children codes are not unique,
	# create unique ids for individuals
##	warning("Sjekk denne!\n")
	new.ids <- f.check.unique.ids( data.in$cov.data )
	id <- new.ids$ids
	pedIndex <- new.ids$pedIndex
	#
	## New, empty list for output data
	gen.data.ut <- vector(mode = "list", length = gen.list.length)
	#
	## Sort families and check coding
	cat( "   Sorting and re-coding families...\n" )

	if( design %in% c( "triad", "cc.triad" ) ){
		# SELECT LINES OF data CORRESPONDING TO EITHER MOTHER, FATHER OR CHILD
		# NOTE THAT DATA LINES NOT CORRESPONDING TO INDIVIDUALS IDENTIFIED IN THE
		# pedIndex FILE WILL NOT BE SELECTED.
		if(any(is.na(id))) stop("There are missing values among the id's", call. = F)
		#
		id.m <- pedIndex[ ,'id.mother' ]
		id.f <- pedIndex[ ,'id.father' ]
		id.c <- pedIndex[ ,'id.child' ]
		if(any(!is.element(c(id.m, id.f, id.c), c(id, NA)))) stop("Something's wrong with the trio id's", call. = F)
		# Track missing
		na.m <- is.na(id.m)
		na.f <- is.na(id.f)
		na.c <- is.na(id.c)
		#
		# Find corresponding positions. NOTE: Some (or all) of the position valuse may be missing, if the corresponding mother, father, or child is missing
		pos.m <- match( id.m, id )
		pos.f <- match( id.f, id )
		pos.c <- match( id.c, id )
		#
		## Loop over input data
		for( i in 1:gen.list.length ){
			cur.gen.data <- data.in$gen.data[[ i ]]
			# 
			## Set up new empty interlaced ff data set
			new.dim <- c( nrow( pedIndex ), 3 * ncol( cur.gen.data ) )
			all.levels <- levels( cur.gen.data[,] )
			cur.gen.data.ut <- ff::ff( dim = new.dim, vmode = .haplinEnv$.vmode.gen.data, levels = all.levels )
			#
			## Fill inn data from original data set
			# Columns containing relevant data in new set
			seq.m <- c( rbind( seq(1, new.dim[2], 6), seq(2, new.dim[2], 6) ) )
			seq.f <- seq.m + 2
			seq.c <- seq.m + 4
			# Insert existing data
			cur.gen.data.ut[ !na.m, seq.m ] <- cur.gen.data[ pos.m[!na.m], ]
			cur.gen.data.ut[ !na.f, seq.f ] <- cur.gen.data[ pos.f[!na.f], ]
			cur.gen.data.ut[ !na.c, seq.c ] <- cur.gen.data[ pos.c[!na.c], ]
			# Insert missing
			cur.gen.data.ut[ na.m, seq.m ] <- NA
			cur.gen.data.ut[ na.f, seq.f ] <- NA
			cur.gen.data.ut[ na.c, seq.c ] <- NA
			#
			## Enter current output data in list
			gen.data.ut[[i]] <- cur.gen.data.ut
		}
		#
		## now - the covariate data
		new.dim <- c( nrow( pedIndex ), 3 * ( ncol( data.in$cov.data ) - 4 ) + 4 )
		cov.data.reordered <- matrix( NA_character_, nrow = new.dim[1], ncol = new.dim[2] )
		chosen.cols <- c( "id.c",colnames( data.in$cov.data )[-(1:4)] )
		cov.data.reordered[,1] <- pedIndex[ ,"family" ]
		cov.data.reordered[, seq(2, new.dim[2], 3)] <- as.matrix( data.in$cov.data[ pos.m,chosen.cols ], mode = "character" )
		cov.data.reordered[, seq(3, new.dim[2], 3)] <- as.matrix( data.in$cov.data[ pos.f,chosen.cols ], mode = "character" )
		cov.data.reordered[, seq(4, new.dim[2], 3)] <- as.matrix( data.in$cov.data[ pos.c,chosen.cols ], mode = "character" )
		
		labs <- c("m", "f", "c")
	} else if( design == "cc" ){
		# SELECT LINES OF data CORRESPONDING TO CHILDREN IN pedIndex FILE
		# NOTE THAT DATA LINES NOT CORRESPONDING TO INDIVIDUALS IDENTIFIED IN THE
		# pedIndex FILE WILL NOT BE SELECTED.
		pos.c <- match( pedIndex[ ,'id.child' ], id )


		for( i in 1:gen.list.length ){
			cur.gen.data <- data.in$gen.data[[ i ]]
			# DIMENSION FOR NEW DATA SET
			new.dim <- c( nrow( pedIndex ), ncol( cur.gen.data ) )
			# NEW INTERLACED DATA SET
			all.levels <- levels( cur.gen.data[,] )

			cur.gen.data.ut <- ff::ff( cur.gen.data[ pos.c, ], dim = new.dim, vmode = .haplinEnv$.vmode.gen.data, levels = all.levels )
			gen.data.ut[[i]] <- cur.gen.data.ut
			# gen.data.ut <- c( gen.data.ut, list( cur.gen.data.ut ) )
		}
		
		chosen.cols <- c( "id.c",colnames( data.in$cov.data )[-(1:4)] )
		cov.data.reordered <- as.matrix( data.in$cov.data[ pos.c,chosen.cols ], mode = "character" )
		cov.data.reordered <- cbind( pedIndex[ ,"family" ], cov.data.reordered )
		
		labs <- "c"
	}
	
	for( i in 1:gen.list.length ){
		cur.gen.data <- data.in$gen.data[[ i ]]
		
		markers1 <- grep( "_a", colnames( cur.gen.data ) )
		markers2 <- grep( "_b", colnames( cur.gen.data ) )

		marker.names.a <- as.vector( t( outer( colnames( cur.gen.data )[ markers1 ], labs, paste, sep = "_" ) ) )
		marker.names.b <- as.vector( t( outer( colnames( cur.gen.data )[ markers2 ], labs, paste, sep = "_" ) ) )
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
