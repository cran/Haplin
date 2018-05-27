f.gwaa.to.ff <- function( data, pedIndex, design = "triad" ){
	## based on gwaaToHaplin, but performing less processing and the output is a list:
	## cov.data - matrix with covariate data
	## gen.data - ff matrix with the genotype data
	## (i.e., the same as output of genDataRead, class "ped.data")
	## this can be read in by genDataPreprocess
	
	n.cols <- dim( data@gtdata )[2]
	n.chunks <- ceiling( n.cols / .haplinEnv$.nb.cols.per.chunk )

	gen.data.list <- list()
	gen.levels <- c()
	for( i in 1:n.chunks ){
		cur.cols <- ( ( i - 1 )*.haplinEnv$.nb.cols.per.chunk + 1 ):min( n.cols, i*.haplinEnv$.nb.cols.per.chunk )
		# this gives the genotype data, alleles separated by "/"
		gen.data.tmp <- as.character( data[ ,cur.cols ] )
		gen.data.split <- f.split.matrix( gen.data.tmp, split = "/" )

		cur.levels <- unique( as.vector( gen.data.split ) )
		recode.na <- FALSE
		if( 0 %in% cur.levels ){
			recode.na <- TRUE
			na.symbol <- 0
		} else if( "0" %in% cur.levels ){
			recode.na <- TRUE
			na.symbol <- "0"
		}
		if( recode.na ){
			gen.data.split[ gen.data.split == na.symbol ] <- NA
			cur.levels[ cur.levels == na.symbol ] <- NA
		}
		gen.levels <- union( gen.levels, cur.levels )

		gen.data.list <- c( gen.data.list, list( ff::as.ff( gen.data.split, vmode = .haplinEnv$.vmode.gen.data, levels = gen.levels ) ) )
	}

	# phenotype data:
	phdata <- phdata( data )
	nph <- ncol( phdata )
	colnames( phdata )[ colnames( phdata ) == "ph" ] <- "cc"
	# RECODE SEX VARIABLE. GenABEL USES male = 1, female = 0, Haplin USES male = 1, female = 2
	phdata$sex <- 2 - phdata$sex

	pedIndex.tab <- read.table( pedIndex, header = T, stringsAsFactors = F )
	cov.colnames <- .haplinEnv$.cov.data.colnames
	phdata <- t( apply( phdata, 1, function(x){
		cur.id <- as.character( x[ "id" ] )
		if( cur.id %in% pedIndex.tab$id.child ){
			idx <- which( pedIndex.tab$id.child == cur.id )
			out.row <- c( pedIndex.tab$family[ idx ], cur.id, pedIndex.tab$id.father[ idx ], pedIndex.tab$id.mother[ idx ], as.character( x[ "sex" ] ), as.character( x[ "cc" ] ) )
		} else {
			idx <- as.numeric( which( pedIndex.tab[ c( "id.mother", "id.father" ) ] == cur.id, arr.ind = TRUE )[ ,1 ] )[ 1 ]
			out.row <- c( pedIndex.tab$family[ idx ], cur.id, "0", "0", as.character( x[ "sex" ] ), as.character( x[ "cc" ] ) )
		}
	} ) )

	colnames( phdata ) <- .haplinEnv$.cov.data.colnames

	out <- list( cov.data = phdata, gen.data = gen.data.list )
	class( out ) <- get( ".class.data.read.ped", envir = .haplinEnv )
	return( out )
}
