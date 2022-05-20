coef.haptable.hapSlide <- function( object, plot.signif.only = FALSE, signif.thresh = 0.05, keep.p.overall = FALSE, poo.p.val = FALSE ){

	coef.list <- lapply( object, coef.haptable, poo.p.val )
	infos <- lapply( coef.list, function( x ){ attr( x, "info" ) } )
	ref.method <- infos[[ 1 ]]$haplos$reference.method
	
	coef.list <- lapply( 1:length( coef.list ), function( x ){
		cur.coef <- coef.list[[ x ]]
		cur.est.names <- rownames( cur.coef )
		rownames( cur.coef ) <- NULL
		tmp.haplos <- object[[ x ]]$haplos
		tmp.haplofreq <- round( object[[ x ]]$haplofreq * 100, digits = 1 )
		if( any( is.na( tmp.haplos ) ) ){
			tmp.haplofreq <- tmp.haplofreq[ !is.na( tmp.haplos ) ]
			tmp.haplos <- tmp.haplos[ !is.na( tmp.haplos ) ]
		}
		cur.haplos <- paste0( tmp.haplos, " (", tmp.haplofreq, "%)" )
		# If the reference method was "ref.cat", remove the reference haplos
		if( ref.method == "ref.cat" ){
			which.ref <- as.numeric( infos[[ x ]]$haplos$ref.cat )
			haplos <- paste0( cur.haplos[ -which.ref ], ",\n ref: ", cur.haplos[ which.ref ] )
			which.ref.rows <- seq( from = which.ref, to = nrow( cur.coef ), by = length( cur.haplos ) )
			cur.coef <- cur.coef[ -which.ref.rows, ]
			cur.haplos <- haplos
			cur.est.names <- cur.est.names[ -which.ref.rows ]
		}
		cbind( marker = names( coef.list )[ x ],
		       haplos = cur.haplos,
		       est.name = cur.est.names,
		       cur.coef,
		       stringsAsFactors = FALSE )
	} )
	# This makes a nice table for plotting with ggplot2
	coefs <- do.call( rbind, coef.list )
	# these give the haplotype frequency, with lower and upper estimate
	# we don't need these, so we remove them from 'coefs'
	which.rows.pvals <- grep( pattern = "^p[[:digit:]]$", as.character( coefs$est.name ) )
	coefs <- coefs[ -which.rows.pvals, ]

	if( keep.p.overall ){
		all.pv.overall <- sapply( object, function(x){
			x$pv.overall[1]
		} )
		all.pv.overall.df <- data.frame(
		  marker = names( all.pv.overall ),
		  haplos = NA,
		  est.name = "pv.overall",
		  est. = all.pv.overall,
		  lower = NA,
		  upper = NA,
		  p.value = all.pv.overall,
		  stringsAsFactors = FALSE )
		coefs <- rbind( coefs, all.pv.overall.df )
	}
	
	# Marking the significant values
	coefs$star <- ""
	coefs$star[ coefs$p.value <= 0.05 ] <- "*"
	coefs$star[ coefs$p.value <= 0.01 ] <- "**"
	coefs$star[ coefs$p.value <= 0.001 ] <- "***"
	
	if( plot.signif.only ){
		all.signif <- coefs$p.value <= signif.thresh
		all.signif.markers <- as.character( unique( coefs$marker[ all.signif ] ) )
		all.matched <- lapply( all.signif.markers, function(x){
			which( coefs$marker == x )
		} )
		all.idx <- Reduce( f = union, x = all.matched )
		coefs.new <- coefs[ all.idx, ]
		if( nrow( coefs.new ) == 0 ){
			warning( "Too low threshold - no data left! Plotting all." )
		} else {
			coefs <- coefs.new
		}
	}
	
	attr( coefs, "info.list" ) <- infos
	return( coefs )
}
