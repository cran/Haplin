#' Pre-processing of the genetic data
#'
#' This function prepares the data to be used in Haplin analysis
#'
#' @param data.in Input data, as loaded by \link{genDataRead} or \link{genDataLoad}.
#' @param map.file Filename (with path if the file is not in current directory) of the
#'   .map file holding the SNP names, if available.
#' @param map.header Logical: does the map.file contain a header in the first row?
#'   Default: FALSE.
#' @param design The design used in the study - choose from:
#'   \itemize{
#'     \item \emph{triad} - (default), data includes genotypes of mother, father and child;
#'     \item \emph{cc} - classical case-control;
#'     \item \emph{cc.triad} - hybrid design: triads with cases and controls
#'   }.
#' @param file.out The core name of the files that will contain the preprocessed data
#'   (character string); ready to load next time with \link{genDataLoad} function;
#'   default: "data_preprocessed".
#' @param dir.out The directory that will contain the saved data; defaults to current 
#'   working directory.
#' @param ncpu The number of CPU cores to use - this speeds up the process for large 
#'   datasets significantly. Default is 1 core, maximum is 1 less than the total number
#'   of cores available on a current machine (even if the number given by the user is
#'   more than that).
#' @param overwrite Whether to overwrite the output files: if NULL (default), will prompt
#'   the user to give answer; set to TRUE, will automatically overwrite any existing files;
#'   and set to FALSE, will stop if the output files exist.
#'
#' @section Details:
#' The .map file should contain at least two columns, where the second one contains SNP 
#'   names. Any additional columns should be separated by a whitespace character, but will 
#'   be ignored. The file should contain a header.
#'
#' @examples
#'   # The argument 'overwrite' is set to TRUE!
#'   # First, read the data:
#'   examples.dir <- system.file( "extdata", package = "Haplin" )
#'   example.file <- file.path( examples.dir, "exmpl_data.ped" )
#'   ped.data.read <- genDataRead( example.file, file.out = "exmpl_ped_data", 
#'    dir.out = tempdir( check = TRUE ), format = "ped", overwrite = TRUE )
#'   ped.data.read
#'   # Take only part of the data (if needed)
#'   ped.data.part <- genDataGetPart( ped.data.read, design = "triad", markers = 10:12,
#'    dir.out = tempdir( check = TRUE ), file.out = "exmpl_ped_data_part", overwrite = TRUE )
#'   # Preprocess as "triad" data:
#'   ped.data.preproc <- genDataPreprocess( ped.data.part, design = "triad",
#'    dir.out = tempdir( check = TRUE ), file.out = "exmpl_data_preproc", overwrite = TRUE )
#'   ped.data.preproc
#'
#' @return  A list object with three elements:
#'   \itemize{
#'     \item \emph{cov.data} - a \code{data.frame} with covariate data (if available in
#'        the input file)
#'     \item \emph{gen.data} - a list with chunks of the genetic data; the data is divided
#'        column-wise, using 10,000 columns per chunk; each element of this list is a
#'        \link[ff]{ff} matrix
#'     \item \emph{aux} - a list with meta-data and important parameters:
#'     \itemize{
#'       \item \emph{variables} - tabulated information of the covariate data;
#'       \item \emph{variables.nas} - how many NA values per each column of covariate data;
#'       \item \emph{alleles} - all the possible alleles in each marker;
#'       \item \emph{alleles.nas} - how many NA values in each marker;
#'       \item \emph{nrows.with.missing} - how many rows contain any missing allele information;
#'       \item \emph{which.rows.with.missing} - vector of indices of rows with missing data (if any)
#'     }.
#'   }
#'
genDataPreprocess <- function(
    data.in = stop( "You have to give the object to preprocess!" ),
    map.file,
    map.header = FALSE,
    design = "triad",
    file.out = "data_preprocessed",
    dir.out = ".",
    ncpu = 1,
    overwrite = NULL
){
	# checking input parameters:
	if( !missing( map.file ) ){
		if( !file.exists( map.file ) ){
			stop( "It seems like the map.file (", map.file, ") doesn't exist! Check and try once more.", call. = FALSE )
		}
	} else if( !is.null( data.in$aux$map.filename ) ){
		map.file <- data.in$aux$map.filename
	} else {
		map.file <- NULL
	}
	
	.info <- data.in$aux$info
	.format <- .info$filespecs$format
	if( !( .format %in% c( "haplin", "ped" ) ) ){
		stop( "Given format (", format,") not recognized! Format has to be one of: ", paste( c("haplin", "ped"), collapse = ", " ), call. = FALSE )
	}
	
	design.list <- get( ".design.list", envir = .haplinEnv )
	if( !( design %in% design.list ) ){
		stop( "Given design(", design,") not recognized! Design has to be one of: ", paste( design.list, collapse = ", " ), call. = FALSE )
	}
	.info$model$design <- design
	
	files.list <- f.make.out.filename( file.out = file.out, dir.out = dir.out, overwrite = overwrite )
	

	if( !is( data.in, "haplin.data" ) ||
	  !all( names( data.in ) == .haplinEnv$.haplin.data.names ) ){
		stop( "The input data is not in the correct format!", call. = FALSE )
	}
	
	if( ncpu < 1 ){
		cat( " You set 'ncpu' to a number less than 1 - resetting it to 1.\n" )
		ncpu <- 1
	}
	#--------------------------------------------

	# creating SNP names (dummy names)
	tot.gen.ncol <- sum( unlist( lapply( data.in$gen.data, ncol ) ) )
	gen.data.colnames <- f.create.snp.names(
	  map.file,
	  header = map.header,
	  ncol = tot.gen.ncol,
	  format = .format,
	  design = design
	)
	marker.names <- gen.data.colnames$marker.names
	gen.data.colnames <- gen.data.colnames$gen.data.colnames

	cur.col <- 1
	for( i in 1:length( data.in$gen.data ) ){
		cur.ncol <- ncol( data.in$gen.data[[ i ]] )
		colnames( data.in$gen.data[[ i ]] ) <- gen.data.colnames[ cur.col:( cur.col + cur.ncol - 1 ) ]
		cur.col <- cur.col + cur.ncol
	}

	# re-organizing data from PED format to internal haplin
	if( .format == "ped" ){
		data.new <- f.ped.to.mfc.new( data.in, design )
	} else {
		data.new <- data.in
		
		if( !is.null( data.new$cov.data ) ){
			orig.colnames <- colnames( data.new$cov.data )
			new.colnames <- paste0( orig.colnames, ".c" )
			colnames( data.new$cov.data ) <- new.colnames
		}
	}
	
	# re-coding the variables per column
	data.recoded <- f.prep.data.parallel( data.new, design, marker.names, ncpu )
	
	## add information
	class( data.recoded ) <- "haplin.ready"
	data.recoded$aux$info <- .info
	data.recoded$aux$class <- class( data.recoded )
	data.recoded$aux$marker.names <- marker.names
	data.recoded$aux$map.filename <- map.file
	
	## find rows with missing data, keep for future reference
	sum.na.list <- lapply( data.recoded$gen.data, function( gen.el.ff ){
		gen.el <- f.extract.ff.numeric( gen.el.ff )
		rowSums( is.na( gen.el ) )
	})
	sum.na <- Reduce( cbind, sum.na.list )
	colnames( sum.na ) <- NULL
	.is.na <- sum.na > 0.1
	rows.with.na <- sum( .is.na )
	data.recoded$aux$nrows.with.missing <- rows.with.na
	data.recoded$aux$which.rows.with.missing <- which( .is.na )
	
	## saving the data in the .RData and .ffData files
	cat( "Saving data... \n" )
	cur.names <- c()
	for( i in 1:length( data.recoded$gen.data ) ){
		cur.name <- paste( get( ".gen.cols.name", envir = .haplinEnv ), i, sep = "." )
		assign( cur.name, data.recoded$gen.data[[i]] )
		cur.names <- c( cur.names, cur.name )
	}
	aux <- data.recoded$aux
	# here, we need to change the format, since now it's 6 columns per 1 marker in the genotype matrix
	aux$info$filespecs$format <- "haplin"
	
	save.list <- c( cur.names, "aux" )
	if( !is.null( data.recoded$cov.data ) ){
		cov.data.in <- data.recoded$cov.data
		save.list <- c( save.list, "cov.data.in" )
	}
	ff::ffsave( list = save.list, file = file.path( dir.out, files.list$file.out.base ) )
	cat( "... saved to files:", files.list$file.out.ff, ", ", files.list$file.out.aux, "\n" )
	return( data.recoded )
}
