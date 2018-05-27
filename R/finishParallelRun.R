#' Closing the Rmpi cluster
#'
#' This function closes all the slaves spawned in \link{initParallelRun} and
#'   finishes the mpi routines. This function MUST BE called after all the \link{haplinSlide}
#'   calls and right before exiting the script/R session!

finishParallelRun <- function(){
	is.spawned <- get( "Rmpi.spawned", envir = .haplinEnv )
	if( !is.spawned ){
		stop( "There are no R slaves to shut down!" )
	}

	Rmpi::mpi.close.Rslaves( dellog = FALSE )
	Rmpi::mpi.finalize( )
	
	assign( "Rmpi.spawned", FALSE, envir = .haplinEnv )
}
