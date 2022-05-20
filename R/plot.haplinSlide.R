#' Plotter function for haplinSlide.
#'
#' This will plot any haplinSlide object in one figure.
#' 
#' The \code{\link{haplinSlide}} object is a list of \code{\link{haplin}} results - by
#'   default in \code{\link{haptable}} form. This is used to plot the relative risk
#'   estimates for all the markers in one plot, to enable easy comparison. NB: those 
#'   estimates that have infinite confidence interval will not be plotted.
#' 
#' @param x The \code{\link{haplinSlide}} object (NB: only the output produced by
#'   running \code{\link{haplinSlide}} with the \code{table.output} argument set to TRUE!)
#' @param filename If the plot should be saved to the disk, give the name of the output
#'   file including the file extension.
#' @param title If the user wishes to override the default title of the plot, give it here.
#' @param windows Numerical vector. If given, this will only plot the chosen windows.
#' @param plot.signif.only Logical: whether to filter out the "non-significant" markers
#'   from the plot. Default: FALSE, i.e., plot everything.
#' @param signif.thresh The threshold defining the significant p-values: if
#'   \code{plot.signif.only == TRUE}, then only the markers with relative risk p-values
#'   lower than the threshold will be kept for plotting. Default: 0.05.
#' @param y.lim Vector with two numbers setting the Y limits of the plotted graph.
#' @param x.title Title for the X axis (default: "marker").
#' @param y.title Title for the Y axis (default: "RR (log scale)").
#' @param file.width Width (in inches) for the output plot, if a filename was given.
#' @param file.height Height (in inches) for the output plot, if a filename was given.
#' @param ... other arguments (ignored).
#' 
#' @return \code{\link[ggplot2]{ggplot}} object.
#'
plot.haplinSlide <- function( x, filename, title, windows, plot.signif.only = FALSE, signif.thresh = 0.05, y.lim, x.title, y.title, file.width, file.height, ... ) {
	
	if( !requireNamespace( "ggplot2", quietly = TRUE ) ){
		stop( "ggplot2 cannot be loaded!", call. = FALSE )
	}
	
	hapSlide.obj <- x
	if( !missing( windows ) ){
		if( length( windows ) > length( x ) | any( !windows %in% 1:length( hapSlide.obj ) ) ){
			stop( "The given window range must be within the haplinSlide result list!", call. = FALSE )
		}
		hapSlide.obj <- x[ windows ]
	}
	
	if( length( hapSlide.obj ) > 50 & !plot.signif.only ){
		stop( "This haplinSlide object is too large - try plotting first only p-values (with 'plotPValues') and deciding which windows you want to choose for plotting with this function.", call. = FALSE )
	}

	info.hapSlide.obj <- attr( hapSlide.obj[[ 1 ]], "info" )
	calc.poo <- info.hapSlide.obj$model$poo
	calc.maternal <- info.hapSlide.obj$model$maternal
		
	# Type of plotting - dependent on the number of haplotypes
	plot.type <- "single"
	if( any( sapply( hapSlide.obj, nrow ) > 2 ) ){
		plot.type <- "multi"
	}

	coefs <- coef.haptable.hapSlide( hapSlide.obj, plot.signif.only, signif.thresh, poo.p.val = calc.poo )
	n.markers.plotted <- length( unique( coefs$marker ) )
	
	# If the PoO effects were calculated, do not plot the separate RR for mother and father,
	# but only the relative RRR_cm_cf
	if( calc.poo ){
		which.rcm <- grep( pattern = "RRcm[[:digit:]]$", x = coefs$est.name )
		which.rcf <- grep( pattern = "^RRcf[[:digit:]]$", x = coefs$est.name )
		coefs <- coefs[ -c( which.rcm, which.rcf ), ]
	}

	# Set the labels, for plotting only
	coefs$labels <- sapply( coefs$est.name, function(x){
		sub( pattern = "RR", replacement = "", x = x )
	} )
	coefs$labels <- sapply( coefs$labels, function(x){
		sub( pattern = "RR", replacement = "", x = x )
	} )
	coefs$labels <- sapply( coefs$labels, function(x){
		sub( pattern = "[[:digit:]]$", replacement = "", x = x )
	} )
	if( calc.maternal ){
		coefs$plot.pos <- sapply( coefs$labels, function(x){
			sub( x = sub( pattern = "^c[[:alnum:]]*$", replacement = "child", x = x ),
					pattern = "^m[[:alnum:]]*$", replacement = "mother" )
		} )
	}

	fig.title <- "Relative risks for child haplotypes"
	if( calc.maternal ){
		fig.title <- "Relative risks for haplotypes"
	}
	if( !missing( title ) ){
		fig.title <- title
	}
	
	# Check if there are very unreliable estimates and remove them - they mess with plotting
	if( any( is.infinite( coefs$lower ) ) | any( is.infinite( coefs$upper ) ) ){
		which.inf <- union( which( is.infinite( coefs$lower ) ), which( is.infinite( coefs$upper ) ) )
		coefs <- coefs[ -which.inf, ]
	}
	if( nrow( coefs ) == 0 ){
		stop( "Too many unreliable estimates - can't plot!", call. = FALSE )
	}

	dodge <- ggplot2::position_dodge( width = 0.6 )
	x.lab <- unique( paste0( coefs$marker, "\n", coefs$haplos ) )
	
	# ggplot creates factors of character columns, ordering them alphabetically
	# but we want a specific order if this is a haplinStrat result:
	#  first, stratum 'all', then strata as numbered
	x.axis.order <- unique(coefs$marker)
	if("all" %in% x.axis.order){
  	coefs$haplos <- factor(
  	    coefs$haplos, levels = unique(coefs$haplos), ordered = TRUE
  	  )
  	coefs$marker <- factor(
  	    coefs$marker, levels = x.axis.order, ordered = TRUE
  	  )
	}
	
	if( missing( y.lim ) ){
		y.lim <- range( c( coefs$upper, coefs$lower ), na.rm = TRUE )
		if( y.lim[ 1 ] < 1/4 ){
			y.lim[ 1 ] <- 1/16
		} else {
			y.lim[ 1 ] <- 1/4
		}
		if( y.lim[ 2 ] > 4 ){
			y.lim[ 2 ] <- y.lim[ 2 ] + 0.5
		} else {
			y.lim[ 2 ] <- 4
		}
	}
	
	if( missing( x.title ) ){
		x.title <- "marker"
	}
	
	if( missing( y.title ) ){
		y.title <- "RR (log scale)"
	}

	if( plot.type == "multi" ){
		basic.plot.multi <- ggplot2::ggplot( coefs, ggplot2::aes_( x = quote( haplos ), y = quote( est. ) ) ) +
			ggplot2::geom_hline( yintercept = 1 ) +
			ggplot2::geom_errorbar( ggplot2::aes_( x = quote( haplos ), ymin = quote( lower ), ymax = quote( upper ), colour = quote( as.factor( labels ) ) ), position = dodge, width = 0.25 ) +
			ggplot2::geom_label( ggplot2::aes( label = labels, fill = as.factor( labels ), fontface = "bold" ), colour = "white", position = dodge ) +
			ggplot2::geom_text( ggplot2::aes_( label = quote( star ), y = quote( upper ), group = quote( est.name ) ), position = dodge, size = 5 ) +
			ggplot2::labs( title = fig.title, y = y.title, x = x.title ) +
			ggplot2::coord_trans( y = "log2", ylim = y.lim ) + 
			ggplot2::theme( legend.position = "none", axis.text.x = ggplot2::element_text( angle = 90 ) )
		if( calc.maternal ){
			final.plot <- basic.plot.multi + ggplot2::facet_grid( plot.pos ~ marker, scales = "free" )
		} else {
			final.plot <- basic.plot.multi + ggplot2::facet_grid( . ~ marker, scales = "free" )
		}
		height.print <- 16
	} else {
		basic.plot.single <- ggplot2::ggplot( coefs, ggplot2::aes_( x = quote( marker ), y = quote( est. ) ) ) +
			ggplot2::geom_hline( yintercept = 1 ) +
			ggplot2::geom_errorbar( ggplot2::aes_( x = quote( marker ), ymin = quote( lower ), ymax = quote( upper ), colour = quote( as.factor( labels ) ) ), position = dodge, width = 0.25 ) +
			ggplot2::geom_label( ggplot2::aes( label = labels, fill = as.factor( labels ), fontface = "bold" ), colour = "white", position = dodge ) +
			ggplot2::geom_text( ggplot2::aes_( label = quote( star ), y = quote( upper ), group = quote( est.name ) ), position = dodge, size = 5 ) +
			ggplot2::labs( title = fig.title, y = y.title, x = x.title ) +
			ggplot2::coord_trans( y = "log2", ylim = y.lim ) + 
			ggplot2::theme( legend.position = "none", axis.text.x = ggplot2::element_text( angle = 90 ) )
		if( calc.maternal ){
			final.plot <- basic.plot.single + ggplot2::scale_x_discrete( labels = x.lab ) + ggplot2::facet_wrap( ~ plot.pos, nrow = 2 )
		} else {
			final.plot <- basic.plot.single + ggplot2::scale_x_discrete( labels = x.lab )
		}
		height.print <- 8
	}
	
	if( !missing( filename ) ){
		ggplot2::ggsave( plot = final.plot, filename = filename, width = min( max( 8, round( n.markers.plotted/2 ) ), 40 ), height = height.print )
	} else {
		print( final.plot )
	}
	
	return( invisible( final.plot ) )
}
