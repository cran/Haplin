f.prep.data.parallel <-  function( data, design, marker.names, ncpu ){
##
## TEST WHETHER FAMILY IS "MFC" OR "C" (USING 6 VERSUS 2 COLUMNS)
if( design %in% c( "triad", "cc.triad" ) ) {
	.t <- 6
}else if( design == "cc" ) .t <- 2

variables <- 0
if( !is.null( data$cov.data ) ){
	variables <- ncol( data$cov.data )
}

if( variables > 0 ){
	cat( "Recoding covariate data...\n" )
	.seq <- 1:variables
	.xdata <- data$cov.data
	.xnamevec <- colnames( data$cov.data )
	
	## GET FREQUENCY COUNT
	.tmp <- lapply( .seq, function(i) f.freq.table( unlist( .xdata[,i] ) ) )
	names( .tmp ) <- .xnamevec
	.freq.x <- lapply( .tmp, function(x) x$tab )
	.unique.x <- lapply( .freq.x, names )
	.nas.x <- sapply( .tmp, function(x) x$nas )
	
	## RECODE VARIABLES DATA, REPLACES ACTUAL CODES WITH 1, 2, 3 ETC, ACCORDING TO .unique.x :
	.xdata.ny <- as.matrix(.xdata)
	.xdata.ny[,] <- NA_integer_
	mode(.xdata.ny) <- "integer"
	#
	for(i in .seq){
		for(j in seq(along = .unique.x[[i]])){
			.xdata.ny[,i][ .xdata[,i] == .unique.x[[i]][j] ] <- j
		}
	}
	.xdata <- as.dframe( .xdata.ny )
	cat( "...done\n" )
}
#
## GENETIC DATA :
nloci <- length( marker.names )
cat( "Recoding genetic data (no. of loci: ", nloci, ")...\n", sep = "" )

## the re-coding of genetic data will be run in parallel
run.para <- FALSE
if( ncpu > 1 ){
	run.para <- TRUE
}
#
## this is the first function that will be called by each worker:
get.levels.snp <- function( my.no, gen.data, tt ){
	# here, 'my.no' is the SNP number
	which.cols <- ( 1 + (my.no - 1)*tt ):( tt*my.no )
	# since the gen.data is in chunks, we need to check in which chunk our data is
	gen.data.tmp <- f.get.gen.data.cols( gen.data, which.cols )

	## FIND THE FREQUENCY OF EACH ALLELE, LIST HAS LENGTH EQUAL TO THE NUMBER OF LOCI:
	.tmp1 <- f.freq.table( gen.data.tmp[,], withNA = TRUE )

	.alleles <- .tmp1$tab
	.nas <- .tmp1$nas

	## RETRIEVE UNIQUE CODES AT EACH MARKER
	.unique <- names( .alleles )
	if( any( is.na( .unique ) ) ){
		.unique <- .unique[ !is.na(.unique) ]
	}
	return( list( list( unique = .unique, allel = .alleles, NAs = .nas ) ) )
}
#
## the recoding function:
recode.snp <- function( my.no, gen.data, unique, tt ){
	# here, my.no is the SNP number
	which.cols <- ( 1 + ( my.no - 1 )*tt ):( tt*my.no )
	gen.list.info <- f.get.which.gen.el( which.cols, ncol( gen.data[[1]] ) )
	gen.data.tmp <- f.get.gen.data.cols( gen.data, which.cols )

	## RECODE GENETIC DATA, REPLACES ACTUAL CODES WITH 1, 2, 3 ETC, ACCORDING TO .unique :
	for( i in 1:length( which.cols ) ){
		col.id <- gen.list.info$col.no[ i ]
		el.id <- gen.list.info$chunk.no[ i ]
		for( j in seq( along = unique[[ my.no ]] ) ){
			tmp.subset <- which(gen.data.tmp[ ,i ] == .subset2( unique, my.no )[ j ])
			data.ny[[ el.id ]][ tmp.subset,col.id ] <- j
		}
	}
	return( TRUE )
}
#
## prepare the "cluster"
if( run.para ){
	if( !requireNamespace( "parallel" ) ){
		stop( "You wanted to run a parallel process without 'parallel' package - something is wrong!" )
	}

	max.ncpu <- parallel::detectCores()
	ncpu <- min( max.ncpu, ncpu )
	cat( "   ...using ", ncpu, " cores\n" )
	cl <- parallel::makePSOCKcluster(getOption("cl.cores", ncpu))
	invisible( parallel::clusterEvalQ( cl, requireNamespace( "ff", quietly = TRUE ) ) )
	invisible( parallel::clusterEvalQ( cl, loadNamespace( "ff" ) ) )
	parallel::clusterExport( cl, c( "f.freq.table", "f.get.which.gen.el", "f.get.gen.data.cols" ), envir = environment() )
	#
	#--------
	# this is just to check whether each node on the newly created "cluster" has access to the data
	open.file.workers <- parallel::clusterEvalQ( cl, names( data ) )
	if( !all( unlist( open.file.workers ) == names( data ) ) ){
		stop( paste( "Problem with opening 'genomic' data file in workers:", which( unlist( open.file.workers ) != names( data$gen.data ) ) ), call. = FALSE )
	}
	#--------
	#
	## we go through all the snps and collect the unique coding per snp
	cat( "   ...checking alleles per SNP...\n" )
	unique.levels.list <- parallel::parSapply( cl, 1:nloci, get.levels.snp, gen.data = data$gen.data, tt = .t )
	names( unique.levels.list ) <- marker.names
	unique.only.list <- lapply( unique.levels.list, function(x) x$unique )
	unique.levels.all <- unique( unlist( unique.only.list ) )
	cat( "   ...done, all alleles:", paste( unique.levels.all, collapse = " " ),"\n" )
	#
	## creating the new data matrix:
	data.ny <- lapply( as.list( 1:length( data$gen.data ) ), function( x ){
		ff::ff( NA, vmode = .haplinEnv$.vmode.gen.data, levels = as.character( c( 1:length( unique.levels.all ), NA ) ), dim = dim( data$gen.data[[x]] ), finalizer = "close" )
	} )
	#
	parallel::clusterExport( cl, c( "data.ny" ), envir = environment() )
	## we go through all the snps again and re-code into the new data matrix:
	cat( "   ...recoding SNPs...\n" )
	recode.results <- parallel::parSapply( cl, 1:nloci, recode.snp, gen.data = data$gen.data, unique = unique.only.list, tt = .t )
	cat( "   ...done\n" )
	parallel::stopCluster( cl )
} else {
	## no parallelization
	cat( "   ...running on only one CPU core! This may take some time...\n" )
	cat( "   ...checking alleles per SNP...\n" )
	unique.levels.list <- sapply( 1:nloci, get.levels.snp, gen.data = data$gen.data, tt = .t )
	names( unique.levels.list ) <- marker.names
	unique.only.list <- lapply( unique.levels.list, function(x) x$unique )
	unique.levels.all <- unique( unlist( unique.only.list ) )
	cat( "   ...done, all alleles:", paste( unique.levels.all, collapse = " " ),"\n" )
	#
	## creating the new data matrix:
	data.ny <- lapply( as.list( 1:length( data$gen.data ) ), function( x ){
		ff::ff( NA, vmode = .haplinEnv$.vmode.gen.data, levels = as.character( c( 1:length( unique.levels.all ), NA ) ), dim = dim( data$gen.data[[x]] ), finalizer = "close" )
	} )
	#
	## we go through all the snps again and re-code into the new data matrix:
	cat( "   ...recoding SNPs...\n" )
	recode.results <- sapply( 1:nloci, recode.snp, gen.data = data$gen.data, unique = unique.only.list, tt = .t )
	cat( "   ...done\n" )
}
#
## remember the colnames in gen.data!
for( i in 1:length( data.ny ) ){
	colnames( data.ny[[ i ]] ) <- colnames( data$gen.data[[ i ]] )
}
#
## creating the output data
if( variables > 0 ){
	.data <- make.ff.data.out( covd = .xdata, gend = data.ny, covd.names = .xnamevec, data.as.is = TRUE )
	.data$aux$variables <- .freq.x
	.data$aux$variables.nas <- .nas.x
} else {
	.data <- make.ff.data.out( gend = data.ny, data.as.is = TRUE )
}

alleles <- lapply( unique.levels.list, function( x ){
	return( x$allel[ !is.na( names( x$allel ) ) ] )
} )
.data$aux$alleles <- alleles
.data$aux$alleles.nas <- lapply( unique.levels.list, function( x ){ x$NAs } )
#
##
return(.data)
}
