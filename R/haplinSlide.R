haplinSlide <- function( data, markers = "ALL", winlength = 1, strata = NULL, table.output = TRUE, cpus = 1, para.env = NULL, slaveOutfile = "", printout = FALSE, verbose = FALSE, ... )
{
##
## RUN HAPLIN ON SLIDING WINDOWS
##
## MERK: BRUKER "#i#" TIL AA KOMMENTERE VEKK DET SOM GIR ADVARSEL MED R CMD check, SIDEN DE RETTE PAKKENE IKKE ER INSTALLERT

## RUNS IN PARALLEL ONLY IF cpus IS NUMERIC > 1
if( !missing( cpus ) ){
	if( !is.numeric( cpus ) ){
		stop( 'The number of CPUs "cpus" must be numeric!', call. = F )
	}
	.run.para <- ( cpus > 1 )
	.para <- get( ".default.para.env", envir = .haplinEnv )
} else {
	.para <- NULL
	.run.para <- F
}

## SELECT PARALLEL MECHANISM ( SOME MAY BE DISABLED )
if( !missing( para.env ) ){
	if( !is.null( para.env  ) ){
		all.para.env <- get( ".para.env", envir = .haplinEnv )
		if( para.env %in% all.para.env ){
			.para <- para.env
		} else {
			cat( "You specified parallel environment that is not recognized:\n" )
			cat( para.env, "\n" )
			cat( "(currently accepted values are:", all.para.env,")\n" )
			.para <- get( ".default.para.env", envir = .haplinEnv )
			cat( "Setting the default environment, '", .para,"' package.\n" )
		}
		# if the Rmpi parallel mechanism was set up, get the CPU number from there
		if( get( "Rmpi.spawned", envir = .haplinEnv ) ){
			# mpi.comm.size returns the number of members in communication channel;
			# this includes the 'master' node, so the number of cpus is 1 less
			cpus <- Rmpi::mpi.comm.size( 1 ) - 1
		}
	}
}

## GET HAPLIN DEFAULTS, OVERRIDE DEFAULTS FOR SPECIFIED ARGUMENTS
cur.call <- match.call( )
.info <- f.catch( cur.call, formals( ) )

design <- .info$model$design

if( winlength > 1 ){
	cat( "\nImportant: Remember that SNPs must be in correct physical ordering\nfor a sliding window approach to be meaningful!\n" )
}

.marker.names <- data$aux$marker.names[ .info$filespecs$markers ]
if( length( .marker.names ) != length( .info$filespecs$markers ) ){
	warning( "Something's wrong with the marker names...", call. = F )
}

## FIND MARKERS INCLUDED IN EACH WINDOW. NB: THEY NOW REFER TO MARKERS IN THE _REDUCED_ FILE, NOT THE ORIGINAL ONE
.slides <- f.windows( markers = .info$filespecs$markers, winlength = winlength )
## USE ORIGINAL MARKER SPECIFICATION AS NAMES
.names <- matrix( data$aux$marker.names[ as.numeric( .slides ) ], ncol = ncol( .slides ) )
.names <- f.create.tag( .names, sep = "-" )

## PREPARE TO RUN ON EACH WINDOW
.nres <- dim( .slides )[1]
orig.args <- f.args.from.info( .info )
attr( data, "sel.markers" ) <- TRUE

## FUNCTION TO RUN HAPLIN ON A SINGLE WINDOW
.f.run <- function( i, args ){
	cat( "Running Haplin on Window '", .names[i], "' (", i, "/", .nres, ")...: ", sep = "" )
	
	args$markers <- .slides[i,]
	args$data <- data
	
	## RUN HAPLIN
	if( is.null( strata ) ){
		.res <- try( do.call( "haplin", args ), silent = T )
		## CHECK IF ERRORS
		if( identical(class( .res ), "try-error") ){
			cat( "RUN FAILED\n" )
		}else{
			cat( "done\n" )
		}
		if(table.output){
			.res <- try(haptable( .res ))
			if( identical(class( .res ), "try-error") ){
				cat( "haptable FAILED\n" )
			}
		}
	}else{
		.res <- try( do.call( "haplinStrat", args ), silent = T )
		## CHECK IF ERRORS
		if( identical(class(.res), "try-error") ){
			cat("RUN FAILED\n")
		}else{
			cat( "done\n" )
		}
		if(table.output){
			.gxe.res <- try(gxe(.res))
			if( identical(class(.gxe.res), "try-error") ){
				cat("RUN FAILED\n")
				.res <- .gxe.res # somewhat ad hoc, but at least picks up the error message
			}else{
				# create haptables and join them
				.res <- try(haptable(.res))
				if( identical(class(.res), "try-error") ){
					cat( "haptable FAILED\n" )
				}else{
					.gxe.res <- try(haptable(.gxe.res))
				 	if( identical(class(.gxe.res), "try-error") ){
						cat( "haptable FAILED\n" )
						.res <- .gxe.res # somewhat ad hoc, but at least picks up the error message
					}else{
						.res <- cbind(.res, .gxe.res)
					}
				}
			}
		}
	}
	return( .res )
}

## DO THE RUN
## EITHER IN SEQUENCE ( SINGLE CPU ) 
## OR IN PARALLEL ( MULTIPLE CPUs )
if( !.run.para ){
	cat( "INFO: No parallel environment selected. Will run in sequential mode.\n" )
	if( !missing( slaveOutfile ) ){
		sink( file = slaveOutfile )
	}
	.res.list <- lapply( seq( length.out = .nres ), .f.run, args = orig.args )
	if( !missing( slaveOutfile ) ) sink( )
	names( .res.list ) <- .names
	
	## CHECK FOR HAPLIN FAILURES
	.errs <- sapply( .res.list, class ) == "try-error"
}else{
	cat( "INFO: Will run in parallel mode using the '", .para,"' package.\n", sep = "" )
	## INITIALIZE CPUS
	if( .para == "parallel" ){
		## parallel
		w <- parallel::makeCluster( cpus, outfile = slaveOutfile )
		invisible( parallel::clusterEvalQ( w, suppressPackageStartupMessages( requireNamespace( "Haplin", quietly = TRUE ) ) ) )
		invisible( parallel::clusterEvalQ( w, suppressPackageStartupMessages( loadNamespace( "Haplin" ) ) ) )
# 		parallel::clusterExport( w, c( ".f.run", "f.sel.markers", "f.get.which.gen.el", "f.get.gen.data.cols", "make.ff.data.out", "f.args.from.info", "f.check.pars", "f.extract.ff.numeric" ), envir = environment( ) )
		on.exit( parallel::stopCluster( w ) )
	} else if( .para == "snow" ){
		## SNOW
		stop( "'snow' package not implemented, please use \"parallel\"!", call. = FALSE )
	} else if( .para == "Rmpi" ){
		is.spawned <- get( "Rmpi.spawned", envir = .haplinEnv )
		if( !is.spawned ){
			stop( "You chose to use 'Rmpi' parallel mechanism but forgot to prepare the cluster. Use 'initParallelRun' function before calling 'haplinSlide' for the first time. After all the calculations are done, right before quitting R, remember to use 'finishParallelRun'!", call. = FALSE )
		}
	
		Rmpi::mpi.bcast.Robj2slave( list( .f.run, f.sel.markers, f.get.which.gen.el, f.get.gen.data.cols, make.ff.data.out, f.args.from.info, f.check.pars, f.extract.ff.numeric ) )
	}
	
	## RUN BY SPLITTING ON CPUS
	if( !missing( slaveOutfile ) ){
		sink( file = slaveOutfile, append = T )
	}
		cat( "\n--- Running haplinSlide using ", cpus, " cpu's ---\n", sep = "" )
	if( !missing( slaveOutfile ) ) sink( )
	
	if( .para == "parallel" ){
		.res.list <- parallel::parLapply( w, seq( length.out = .nres ), .f.run, args = orig.args )
	} else if( .para == "Rmpi" ){
		.res.list <- Rmpi::mpi.parLapply( seq( length.out = .nres ), .f.run, args = orig.args, job.num = .nres )
	}
	names( .res.list ) <- .names
	
	## CHECK FOR HAPLIN FAILURES
	.ftest <- function( x ){
		return( identical( class( x ), "try-error" ) )
	}
	if( .para == "parallel" ){
		.errs <- unlist( parallel::parLapply( w, .res.list, .ftest ) )
	} else if( .para == "Rmpi" ){
		.errs <- unlist( Rmpi::mpi.parLapply( .res.list, .ftest ) )
	}
}

## GO THROUGH POSSIBLE ERRORS
if( any( .errs ) ){
	## COLLECT ERROR MESSAGES, CHANGE OUTPUT TO NA WITH ERR. MESS. AS ATTRIBUTES
	.mess <- .res.list[.errs]
	.err.res <- lapply( .mess, function( x ) {
		.tmpres <- NA
		attributes( x ) <- NULL # REMOVE "try-error"
		attr( .tmpres, "error.message" ) <- x
		return( .tmpres )
	} )
	.res.list[.errs] <- .err.res
}
if( !missing( slaveOutfile ) ){
	sink( file = slaveOutfile, append = T )
}
cat( "\n --- haplinSlide has completed ---\n" )
if( !missing( slaveOutfile ) ) sink( )

## SET CLASS
class( .res.list ) <- "haplinSlide"

## RETURN RESULT LIST
return( .res.list )
}
