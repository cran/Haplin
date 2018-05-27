#' Plotter function for haplinStrat results.
#'
#' This will automatically plot all haplinStrat results on one figure.
#'
#' This function uses the same style as \code{\link{plot.haplinSlide}} and plots all of
#'  the \code{\link{haplinStrat}} results on one figure, for easy comparison. NB: those 
#'   estimates that have infinite confidence interval will not be plotted.
#'
#' @inheritParams plot.haplinSlide
#' @param x The \code{\link{haplinSlide}} object (NB: only the output produced by
#'   running \code{\link{haplinSlide}} with the \code{table.output} argument set to TRUE!)
#' @param y.lim Vector with two numbers setting the Y limits of the plotted graph.
#' @param x.title Title for the X axis (default: "marker").
#' @param y.title Title for the Y axis (default: "RR (log scale)").
#' @param file.width Width (in inches) for the output plot, if a filename was given.
#' @param file.height Height (in inches) for the output plot, if a filename was given.
#' @param ... other arguments (ignored).
#'
#' @return \code{\link[ggplot2]{ggplot}} object.
#'
plot.haplinStrat <- function( x, filename, title, windows, plot.signif.only = FALSE, signif.thresh = 0.05, y.lim, x.title, y.title, file.width, file.height, ... ) {
	if( !requireNamespace( "ggplot2", quietly = TRUE ) ){
		stop( "ggplot2 cannot be loaded!", call. = FALSE )
	}
	
	plot.obj <- lapply( x, haptable )
	if( missing( x.title ) ){
		x.title <- "stratum"
	}
	final.plot <- plot.haplinSlide( plot.obj, filename, title, windows, plot.signif.only, signif.thresh, y.lim, x.title, y.title, file.width, file.height )
	
	return( invisible( final.plot ) )
}
