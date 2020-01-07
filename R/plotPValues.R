#' Plotting p-values for relative risks
#'
#' This function plots p-values for the relative risks calculated by \code{\link{haplinSlide}}.
#' 
#' The output of \code{\link{haplinSlide}} can be very lengthy and not suitable for
#'   an overall plot of all the relative risks (RR) on one figure. Therefore, it's advised
#'   to first plot only the p-values for each window (user can choose which p-values to
#'   plot - see parameter \code{which.p.val}), and only then plot the RRs for specific
#'   windows, for which the p-values are significant.
#'
#' @param object The \code{\link{haplinSlide}} results: list of \code{\link{haptable}}
#'   objects.
#' @param windows Numerical vector; if given, the plot will be restricted to only those.
#' @param which.p.val Character string specifying which p-values to choose for plotting:
#'   "overall" (default), "child", "child.double", "maternal", "maternal.double",
#'   "paternal". The last three options can be chosen only if \code{\link{haplinSlide}}
#'    was run with \code{maternal = TRUE} or \code{poo = TRUE}.
#' @param plot.signif.only Logical: whether to filter out the "non-significant" markers
#'   from the plot. Default: FALSE, i.e., plot everything.
#' @param signif.thresh The threshold defining the significant p-values: if
#'   \code{plot.signif.only == TRUE}, then only the markers with relative risk p-values
#'   lower than the threshold will be kept for plotting. Default: 0.05.
#' @param title Optional character string for the title of the figure.
#' @param filename If the plot should be saved to the disk, give the name of the output
#'   file including the file extension.
#'
#' @return Invisibly returns the table with only the plotted p-values.
#'
plotPValues <- function( object, windows, which.p.val = "overall", plot.signif.only = FALSE, signif.thresh = 0.05, title, filename ) {
	if( !requireNamespace( "ggplot2", quietly = TRUE ) ){
		stop( "ggplot2 cannot be loaded!", call. = FALSE )
	}
	
	hapSlide.obj <- object
	if( !missing( windows ) ){
		if( length( windows ) > length( object ) | any( !windows %in% 1:length( hapSlide.obj ) ) ){
			stop( "The given window range must be within the haplinSlide result list!", call. = FALSE )
		}
		hapSlide.obj <- object[ windows ]
	}
	
	all.p.val.names <- c( "overall", "child", "child.double", "maternal", "maternal.double", "paternal", "poo" )
	if( !which.p.val %in% all.p.val.names ){
		stop( "'which.p.val' should be one of: ", paste( all.p.val.names, collapse = "," ), call. = FALSE )
	}
	
	# ---
	# construct a "map" to find the correct p-values in coefs table
	all.coefs.p.val <- vector( mode = "list", length = length( all.p.val.names ) )
	names( all.coefs.p.val ) <- all.p.val.names
	all.coefs.p.val$overall <- "pv.overall"
	all.coefs.p.val$child <- c( "RRc", "RR", "RR.est." )
	all.coefs.p.val$child.double <- c( "RRdd", "RRcdd" )
	all.coefs.p.val$maternal <- c( "RRm", "RRcm" )
	all.coefs.p.val$maternal.double <- "RRmdd"
	all.coefs.p.val$paternal <- c( "RRf", "RRcf" )
	all.coefs.p.val$poo <- "RRcm_RRcf"
	# ---

	# need to remove the runs that resulted in error:
	which.obj.na <- which( sapply( hapSlide.obj, function(x){ !is.data.frame( x ) } ) )
	if( length( which.obj.na ) != 0 ){
		hapSlide.obj <- hapSlide.obj[ -which.obj.na ]
	}
	
	poo.p.val <- FALSE
	if( which.p.val == "poo" ){
		poo.p.val <- TRUE
	}
	coefs <- coef.haptable.hapSlide( hapSlide.obj, plot.signif.only, signif.thresh, keep.p.overall = TRUE, poo.p.val = poo.p.val )
	
	unique.coefs.p.val <- unique( as.character( coefs$est.name ) )
	chosen.p.val.names <- all.coefs.p.val[[ which.p.val ]]
	if( which.p.val == "overall" ){
		name.patterns <- "pv.overall"
	} else {
		name.patterns <- paste0( "^", chosen.p.val.names, "[[:digit:]*]" )
	}
	for( pattern in name.patterns ){
		which.p.val.plot <- grep( pattern = pattern, x = as.character( coefs$est.name ) )
		if( length( which.p.val.plot ) != 0 ){
			break
		}
	}
	if( length( which.p.val.plot ) == 0 ){
		stop( "There are no ", which.p.val, " estimates! Check the data and try once more.", call. = FALSE )
	}
	
	coefs.plotting <- coefs[ which.p.val.plot, ]
	if( !all( is.na( coefs.plotting$haplos ) ) ){
		coefs.plotting$haplos <- paste0( coefs.plotting$marker, ":\n", coefs.plotting$haplos )
	} else {
		coefs.plotting$haplos <- coefs.plotting$marker
	}
	
	if( plot.signif.only ){
		# remove the non-significant ones:
		which.signif <- coefs.plotting$p.value <= signif.thresh
		coefs.plotting <- coefs.plotting[ which.signif, ]
	}
	
	if( nrow( coefs.plotting ) == 0 ){
		stop( "There are no ", which.p.val, " estimates lower than ", signif.thresh, "!", call. = FALSE )
	}

	fig.title <- paste0( which.p.val, " p-values" )
	if( !missing( title ) ){
		fig.title <- title
	}

	p.val.plot <- ggplot2::ggplot( coefs.plotting, ggplot2::aes_( quote( haplos ), quote( p.value ) ) ) +
		ggplot2::geom_point( ggplot2::aes_( colour = quote( p.value ) ) ) +
		ggplot2::geom_hline( yintercept = signif.thresh, linetype = 2 ) +
		ggplot2::scale_color_gradient( low = "blue" ) +
		ggplot2::theme( axis.text.x = ggplot2::element_text( angle = 90 ), axis.title.x = ggplot2::element_blank() )
	
	n.markers.plotted <- nrow( coefs.plotting )
	
	if( !missing( filename ) ){
		ggplot2::ggsave( filename = filename, width = min( max( 8, round( n.markers.plotted/2 ) ), 40 ), height = 8 )
	} else {
		print( p.val.plot )
	}
	
	return( invisible( coefs.plotting ) )
}















