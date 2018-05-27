#' Extracting part of genetic data.
#'
#' This function enables to extract (and save for later use) part of genetic data
#' read in with \link{genDataRead}.
#'
#' The genetic data from GWAS studies can be quite large, and thus the analysis
#' is time-consuming. If a user knows where they want to focus the analysis,
#' they can use this function to extract part of the entire dataset and use
#' only this part in subsequent Haplin analysis.
#'
#' @param data.in The data object (in format as the output of \link{genDataRead}).
#' @param file.out The base for the output filename (default: "my_data_part").
#' @param dir.out The path to the directory where the output files will be saved.
#' @param design The design used in the study - choose from:
#'   \itemize{
#'     \item \emph{triad} - (default), data includes genotypes of mother, father and child;
#'     \item \emph{cc} - classical case-control;
#'     \item \emph{cc.triad} - hybrid design: triads with cases and controls;
#'   }.
#'
#' Any of the following can be given to narrow down the dataset:
#' @param markers Numeric vector with numbers indicating which markers to choose.
#' @param indiv.ids Character vector giving IDs of individuals. \strong{CAUTION:}
#'   in a standard PED file, individual IDs are not unique, so this will select
#'   all individuals with given IDs.
#' @param rows Numeric vector giving the positions - this will select only these rows.
#' @param cc One or more values to choose based on case-control status ('cc' column).
#' @param sex One or more values to choose based on the 'sex' column.
#' @param overwrite Whether to overwrite the output files: if NULL (default), will prompt
#'   the user to give answer; set to TRUE, will automatically overwrite any existing files;
#'   and set to FALSE, will stop if the output files exist.
#' @param ... If any additional covariate data are available in \code{data.in}, 
#'   the user can choose based on values of these (see the Examples section).
#'
#' @return A list object with three elements:
#'   \itemize{
#'     \item \emph{cov.data} - a \code{data.frame} with covariate data (if available in
#'        the input file)
#'     \item \emph{gen.data} - a list with chunks of the genetic data; the data is divided
#'        column-wise, using 10,000 columns per chunk; each element of this list is a
#'        \link[ff]{ff} matrix
#'     \item \emph{aux} - a list with meta-data and important parameters.
#'   }
#'   This now contains only the selected subset of data.
#'
#' @examples
#'   # The argument 'overwrite' is set to TRUE!
#'   # Read the data:
#'   examples.dir <- system.file( "extdata", package = "Haplin" )
#'   example.file <- paste0( examples.dir, "/HAPLIN.trialdata2.txt" )
#'   my.gen.data.read <- genDataRead( file.in = example.file, file.out = "trial_data",
#'    dir.out = ".", format = "haplin", allele.sep = "", n.vars = 2, cov.header = 
#'    c( "smoking", "sex" ), overwrite = TRUE )
#'   my.gen.data.read
#'   # Extract part with only men:
#'   men.subset <- genDataGetPart( my.gen.data.read, design = "triad", sex = 1,
#'     dir.out = ".", file.out = "gen_data_men_only", overwrite = TRUE )
#'   men.subset
#'   # Extract the part with only smoking women:
#'   women.smoke.subset <- genDataGetPart( my.gen.data.read, design = "triad",
#'     sex = 0, smoking = c( 1,2 ), overwrite = TRUE )
#'   women.smoke.subset 
#'
#' @section Warning:
#' No checks are performed when choosing a subset of the data - it is the user's
#' obligation to check whether the data subset contains correct number of individuals
#' (especially important when using the \code{triad} design study) and/or markers!
#'

genDataGetPart <- function( data.in = stop( "No data given!", call. = FALSE ), design = stop( "Design type must be given!" ), markers, indiv.ids, rows, cc, sex, file.out = "my_data_part", dir.out = ".", overwrite = NULL, ... ){
	#---- checking the input params ---------------------
	files.list <- f.make.out.filename( file.out = file.out, dir.out = dir.out, overwrite = overwrite )
	
	if( (class(data.in) != "haplin.data") || !all( names( data.in ) == .haplinEnv$.haplin.data.names ) ){
		stop( "The input data is not in the correct format!", call. = FALSE )
	}
	
	design.list <- get( ".design.list", envir = .haplinEnv )
	if( !( design %in% design.list ) ){
		stop( "Given design(", design,") not recognized! Design has to be one of: ", paste( design.list, collapse = ", " ), call. = FALSE )
	}
	#---- done checking the input params ----------------
	#----------------------------------------------------
	# get all arguments
	all.args <- match.call()[-1]
	all.args.names <- names( all.args )
	# arguments that do not define selection of data subset:
	non.sel.args <- c( "data.in", "file.out", "dir.out", "overwrite", "design" )
	
	which.non.sel.args.given <- c()
	for( i in 1:length( non.sel.args ) ){
		which.cur.arg <- match( non.sel.args[ i ], all.args.names )
		if( !is.na( which.cur.arg ) ){
			which.non.sel.args.given <- c( which.non.sel.args.given, which.cur.arg )
		}
	}
	# arguments that define selection of data subset:
	selection.args <- all.args[ -which.non.sel.args.given ]

	# this is the maximum set of arguments for subset selection, based on the read in covariate data
	format <- data.in$aux$info$filespecs$format
	if( format == "ped" ){
		correct.sel.args <- union( names( formals() ), colnames( data.in$cov.data )[ -(1:4) ] )
	} else { # haplin file
		correct.sel.args <- union( names( formals() ), colnames( data.in$cov.data ) )
	}
	correct.sel.args <- correct.sel.args[ -( match( c( non.sel.args, "..." ),
		correct.sel.args ) ) ]
	
	max.rows <- nrow( data.in$gen.data[[1]] )
	all.rows <- 1:max.rows
	all.markers <- sum( unlist( lapply( data.in$gen.data, ncol ) ) )
	all.cols <- 1:all.markers
	
	subset.rows <- all.rows
	subset.cols <- all.cols
	is.subset.cols <- FALSE
	is.subset.rows <- FALSE
	family <- "c"
	if( design %in% c( "triad", "cc.triad" ) & format == "haplin" ){
		family <- "mfc"
	}
	
	cat( "Provided arguments:\n" )
	for( arg in names( selection.args ) ){
		if( arg == "markers" ){
			val <- eval( all.args$markers )
			print.val <- paste( val, collapse = ", " )
			if( length( val ) > 20 ){
				print.val <- paste( paste( head( val ), collapse = ", "), "...", paste( tail( val ), collapse = ", " ) )
			}
			cat( " --- chosen markers:", print.val, "\n" )
			if( any( val > all.markers ) | any( val < 0 ) ){
				stop( "wrong markers chosen; not in range of given data: max ", all.markers, "!\n", call. = FALSE )
			}
			if( identical( val, all.cols ) ){
				warning( "this selection is equal to choosing all markers!", call. = FALSE )
			} else {
				subset.cols <- f.sel.markers( n.vars = 0, markers = val, family = family, split = TRUE, ncols = all.markers )
				is.subset.cols <- TRUE
			}
			
		} else if( arg == "indiv.ids" ){
			val <- eval( all.args$indiv.ids )
			print.val <- paste( val, collapse = ", " )
			if( length( val ) > 20 ){
				print.val <- paste( paste( head( val ), collapse = ", "), "...", paste( tail( val ), collapse = ", " ) )
			}
			cat( " --- individual IDs:", print.val, "\n" )
			
			if( is.null( data.in$cov.data ) |
				!( "id.c" %in% colnames( data.in$cov.data ) ) ){
				# this can happen when the haplin-formatted file had no extra info
				stop( "There is no information on individual IDs in the data!", call. = FALSE )
			}
			
			chosen.indiv.logic <- val %in% as.character( data.in$cov.data[ ,"id.c" ] )
			if( !any( chosen.indiv.logic ) ){
				stop( "wrong individual IDs chosen!", call. = FALSE )
			} else if( sum( chosen.indiv.logic ) == max.rows ){
				warning( "this selection is equal to choosing all the individuals available in dataset!", call. = FALSE )
			} else {
				chosen.indiv <- do.call( c, sapply( val, function( x ){
					list( which( as.character( data.in$cov.data[ ,"id.c" ] ) == x ) )
				} ) )
				subset.rows <- intersect( subset.rows, chosen.indiv[ !( is.na( chosen.indiv ) ) ] )
				is.subset.rows <- TRUE
			}
			
		} else if( arg == "rows" ){
			val <- eval( all.args$rows )
			print.val <- paste( val, collapse = ", " )
			if( length( val ) > 20 ){
				print.val <- paste( paste( head( val ), collapse = ", "), "...", paste( tail( val ), collapse = ", " ) )
			}
			cat( " --- rows chosen:", print.val, "\n" )
			if( any( val < 1 ) | any( val > nrow( data.in$cov.data ) ) ){
				stop( "wrong rows, not in range of given data: max ", max.rows, "\n", call. = FALSE )
			}
			if( identical( val, all.rows ) ){
				warning( "this selection is equal to choosing all the rows available in dataset!", call. = FALSE )
			} else {
				subset.rows <- intersect( subset.rows, val )
				is.subset.rows <- TRUE
			}
			
		} else if( !( arg %in% correct.sel.args ) ){
 			stop( "the given argument: ", arg, " is not recognizable!", call. = FALSE )
 			
		} else { #any other argument if additional covariate data are loaded
			val <- eval( all.args[[ match( arg, names( all.args ) ) ]] )
			print.val <- paste( val, collapse = ", " )
			if( length( val ) > 20 ){
				print.val <- paste( paste( head( val ), collapse = ", "), "...", paste( tail( val ), collapse = ", " ) )
			}
			cat( " --- argument: ", arg,", chosen values: ", print.val, "\n", sep = "" )
			avail.vals <- unique( as.numeric( data.in$cov.data[ ,arg ] ) )
			if( !all( val %in% avail.vals ) ){
				stop( "wrong argument(", arg, ") chosen, available: ", paste( avail.vals, sep = " ", collapse = " " ), "!", call. = FALSE )
			}
			if( all( avail.vals %in% val ) ){
				warning( "this selection is equal to choosing all the possible values of the given covariate (", arg, ")!", call. = FALSE )
			} else {
				chosen.indiv <- do.call( c, sapply( val, function( x ){
					list( which( as.numeric( data.in$cov.data[ ,arg ] ) == x ) )
				} ) )
				subset.rows <- intersect( subset.rows, chosen.indiv )
				is.subset.rows <- TRUE
			}
			
	 	}
	}
	if( is.subset.rows | is.subset.cols ){
		cat( "INFO: Will select ", length( subset.rows ), " rows and ", length( subset.cols ), " columns.\n", sep = "" )
	} else {
		stop( "No subset selected.\n", call. = FALSE )
	}
	
	if( is.subset.cols ){
		# check how many chunks will be needed
		nb.cols.per.chunk <- get( ".nb.cols.per.chunk", envir = .haplinEnv )
		nb.chunks <- ceiling( length( subset.cols ) / nb.cols.per.chunk )
		
		gen.data.col.wise <- sapply( 1:nb.chunks, function( x ){
			first.col <- ( nb.cols.per.chunk*( x - 1 ) + 1 )
			last.col <- min( ( nb.cols.per.chunk*x ), length( subset.cols ) )
			cur.cols <- subset.cols[ first.col:last.col ]
			list( f.get.gen.data.cols( data.in$gen.data, cur.cols ) )
		} )
	} else {
		gen.data.col.wise <- data.in$gen.data
	}
	
	cov.data.in <- NULL
	if( is.subset.rows ){
		# need to choose from both gen.data and cov.data
		gen.data.col.wise <- lapply( gen.data.col.wise, function( x ){
			sub <- x[ subset.rows, ]
			ff::ff( sub, levels = ff::levels.ff( sub ), dim = dim( sub ), vmode = ff::vmode( x ) )
		})
		
		if( !is.null( data.in$cov.data ) ){
			cov.data.in <- data.in$cov.data[ subset.rows, ]
		}
	} else if( !is.null( data.in$cov.data ) ){
		cov.data.in <- data.in$cov.data
	}
	data.out <- list( cov.data = cov.data.in, gen.data = gen.data.col.wise, aux = data.in$aux )
	class( data.out ) <- class( data.in )
	
	## saving the chosen part of the data
	cat( "Saving data... \n" )
	cur.names <- c()
	for( i in 1:length( gen.data.col.wise ) ){
		cur.name <- paste( get( ".gen.cols.name", envir = .haplinEnv ), i, sep = "." )
		assign( cur.name, gen.data.col.wise[[i]] )
		cur.names <- c( cur.names, cur.name )
	}
	aux <- data.in$aux
	save.list <- c( cur.names, "aux" )
	if( !is.null( cov.data.in ) ){
		save.list <- c( save.list, "cov.data.in" )
	}
	ff::ffsave( list = save.list, file = paste0( dir.out, "/", files.list$file.out.base ) )
	cat( "... saved to files: ", files.list$file.out.ff, ", ", files.list$file.out.aux, "\n", sep = "" )

	return( data.out )
}
