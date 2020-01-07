#' Initialization of the Rmpi cluster
#'
#' This function prepares a cluster using \code{Rmpi} package. The initialization
#'   is paired with closing of the cluster using the \code{\link{finishParallelRun}} function.
#'
#' @param cpus Number of cores to use (optional). If given, only that number of CPUs 
#'   will be used. By default (if not set), the Rmpi will check how many CPUs are
#'   available in the system and take the maximum number.

initParallelRun <- function( cpus ){

	if( !requireNamespace( "Rmpi", quietly = TRUE ) ){
		stop( "Rmpi cannot be loaded!" )
	}
	
	if( Rmpi::mpi.universe.size( ) > 1 ){
#	cat( "MPI universe size: ", mpi.universe.size(), "\n" )
		# spawn as many as possible
		cpus.max <- Rmpi::mpi.universe.size() - 1
	} else {
		cpus.max <- parallel::detectCores()
	}
	if( cpus.max < 2 ){
		stop( "There is only 1 core available. No need to run parallel!" )
		assign( "Rmpi.spawned", FALSE, envir = .haplinEnv )
	}

	if( missing( cpus ) ){
		cpus <- cpus.max
	} else {
		if( cpus > cpus.max ){
			cat( "You requested ", cpus, " CPUs, while your system seems to have only ", cpus.max, " available. Correcting.\n" )
			cpus <- cpus.max
		}
	}
	
	Rmpi::mpi.spawn.Rslaves( nslaves = cpus )
	
	cat( "\nINFO: Running in 'Rmpi' parallel mode, on ", cpus, "CPUs.\n\n" )
	
	# if R closes unexpectedly, cleanup!
	.Last <- function( ){
		if( is.loaded( "mpi_initialize" ) ){
			if( Rmpi::mpi.comm.size( 1 ) > 0 ){
				Rmpi::mpi.close.Rslaves( dellog = FALSE )
			}
			Rmpi::mpi.finalize( )
		}
	}
	
	Rmpi::mpi.bcast.cmd( suppressPackageStartupMessages( requireNamespace( "Haplin", quietly = TRUE ) ) )
	Rmpi::mpi.bcast.cmd( suppressPackageStartupMessages( loadNamespace( "Haplin" ) ) )
	Rmpi::mpi.bcast.Robj2slave( list( f.sel.markers, f.get.which.gen.el, f.get.gen.data.cols, make.ff.data.out, f.args.from.info, f.check.pars ) )
	
	assign( "Rmpi.spawned", TRUE, envir = .haplinEnv )
}
