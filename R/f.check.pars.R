f.check.pars <- function( mcall, defaults ){
##
## SET DEFAULTS AND CHECK THE LOGIC OF THE INPUT PARAMETERS TO HAPLIN
## NOTE: DOES NOT CHECK FOR SUPERFLUOUS ARGUMENTS
##
#
## REPLACE DEFAULTS WITH SPECIFIED PARAMETERS
params <- defaults
params[names( mcall )] <- mcall
#

## COLLECT IN LIST
.info <- list( 
	filename = params[["filename"]],
	data = list(),
	filespecs = params[c( "markers", "n.vars", "sep", "allele.sep", "na.strings" )], 
	model = params[c( "design", "use.missing", "xchrom", "comb.sex", "maternal", "poo", "test.maternal", "scoretest" )], 
	variables = params[c( "ccvar", "covar", "strata", "sex" )],
	haplos = params[c( "reference", "response", "threshold", "max.haplos", "haplo.file" )],
	control = params[c( "resampling", "max.EM.iter", "data.out", "verbose", "printout", "sel.markers" )]
	 )
class( .info ) <- "info"
#
## DEFINE SPECIFIC VARIABLES
.xchrom <- .info$model$xchrom
.ccdesign <- .info$model$design %in% c( "cc.triad", "cc" )
#
## BASIC RECURSIVE CHECK FOR CORRECT ARGUMENT VALUES. 
## NOTE: DOES NOT CHECK NULL VALUES
.allowed <- .info
.allowed$model$scoretest <- c( "yes", "no", "only" )
.allowed$model$design <- c( "triad", "cc.triad", "cc" )
.allowed$haplos$response <- c( "mult", "free" )
.allowed$control$data.out <- c( "no", "basic", "prelim", "null", "full" )
.allowed$control$resampling <- c( "no", "jackknife" )
.allowed$model$comb.sex <- c( "males", "females", "single", "double" )
#
for( i in seq( along = .info )[-1] ){# NO CHECK FOR filename
	for( j in seq( along = .info[[i]] ) ){
		if( length( .info[[c( i,j )]] ) == 1 ){# EXCLUDES NULL VALUES AND VECTOR VALUES ( LIKE markers )
			if( !is.element( .info[[c( i,j )]], .allowed[[c( i,j )]] ) ){
				stop( paste( "The argument ", names( .info[[i]] )[j], " has an invalid value. \nIt should be one of: \n", paste( .allowed[[c( i, j )]], collapse = ", " ), sep = "" ), call. = F )
			}
		}
	}
}
#
##
if( .info$model$design %in% c( "triad", "cc.triad" ) ) {
	.info$model$fam <- "mfc"
}
if( .info$model$design == "cc" ){
	.info$model$fam <- "c"
}
#

## SET VALUES
# REPLACE "ALL" IF NECESSARY
.info$control$sel.markers <- TRUE
if( identical( .info$filespecs$markers, "ALL" ) ){
	if( .info$model$design %in% c( "triad", "cc.triad" ) ){
		ncol.per.marker <- 6
	} else if( .info$model$design == "cc" ) {
		ncol.per.marker <- 2
	} else {
		stop( "Given design: ", .info$model$design, " is not recognized!", call. = FALSE )
	}
	
	if( class( mcall$data ) == "haplin.ready" ){
		.info$filespecs$markers <- seq( length.out = sum( unlist( lapply( mcall$data$gen.data, ncol ) ) )/ncol.per.marker )
		attr( mcall$data, "sel.markers" ) <- FALSE
	}
	.info$control$sel.markers <- FALSE
}

# SET CORRECT VALUES FOR COLUMN NUMBERS
.info$filespecs$n.vars <- 0
if( class( mcall$data ) == "haplin.ready" ){
	if( !is.null( mcall$data$cov.data ) ){
		.info$filespecs$n.vars <- ncol( mcall$data$cov.data )
		if( .xchrom ){
			if( any( .haplinEnv$sex.column %in% colnames( mcall$data$cov.data ) ) ){
				.info$variables$sex <- match( .haplinEnv$sex.column, colnames( mcall$data$cov.data ) )
			} else {
				.info$variables$sex <- mcall$sex
# 				stop( "The 'xchrom' parameter is set but the data does not contain information about sex!", call. = FALSE )
			}
		}
		if( .ccdesign ){
			# Get ccvar from haplin argument if present, 
			# else use cc.c from cov.data
			if(is.numeric(mcall$ccvar)){
				.ccvar <- mcall$ccvar
			}else if(is.character(mcall$ccvar)){
				.ccvar <- match( mcall$ccvar, colnames( mcall$data$cov.data ) )
			} else {
				.ccvar <- match( .haplinEnv$cc.column, colnames( mcall$data$cov.data ) )
			}
			if(is.na(.ccvar) || (.ccvar <= 0.01) || (.ccvar >= ncol(mcall$data$cov.data) + 0.01)) stop("'ccvar' must specify a column name or position in the covariate data frame", call. = FALSE)
			.info$variables$ccvar <- .ccvar
		}
	}
}

## ONLY USE response = "mult", AND DISALLOW poo, FOR MALES ONLY IN xchrom
if( .xchrom && !is.null( .info$model$comb.sex ) && .info$model$comb.sex == "males" ){
	if( .info$haplos$response != "mult" ){
		warning( 'Can only use response = "mult" with comb.sex = "males". Has been changed to "mult".', call. = F )
		.info$haplos$response <- "mult"
	}
	if( .info$model$poo ){
		stop( 'parent-of-origin estimation not possible when comb.sex = "males".', call. = F )
	}
}
## FURTHER RESTRICT USE OF PARENT-OF-ORIGIN
if( .info$model$poo ){
	if( .info$model$design == "cc" ){
		stop( 'parent-of-origin effects not available when design = "cc"', call. = F )
	}
	if( .info$haplos$reference == "reciprocal" ){
		warning( 'Can only (for the time being) use reference = "ref.cat" or "population" when poo == TRUE. Has been changed to "ref.cat".', call. = F )
		.info$haplos$reference <- "ref.cat"
	}
}
if((.info$model$design == "cc") & .xchrom) stop('The X-chromosome analysis with a "pure" cc design needs\n some additional testing before being made available again...')
#
## TEST FOR CC VARIABLE SPECIFICATION
if( .ccdesign ){
	if( is.null( .info$variables$ccvar ) ){
		stop( 'Parameter "ccvar" must be specified when using design "cc.triad" or "cc"!', call. = F )
	}
}else{
	if( !is.null( .info$variables$ccvar ) )stop( 'Parameter "ccvar" should only be specified when using design "cc.triad" or "cc"!', call. = F )
}
#
## MAKE SURE sex COLUMN IS SPECIFIED IF ON THE xchrom
if( !is.logical( .xchrom ) ){
	stop( 'Argument "xchrom" must be a logical ( either "TRUE" or "FALSE" )', call. = F )
}
if( F & .xchrom ){
	if( .info$filespecs$n.vars == 0 ){
		stop( 'Argument "n.vars" must be at least 1 to allow for a sex variable when "xchrom = TRUE"', call. = F )
	}
	if( !is.numeric( .info$variables$sex ) ){
		stop( 'Argument "sex" should be a numeric value ( the column number of the sex variable ) when "xchrom = TRUE"', call. = F )
	}
	if( .info$variables$sex > .info$filespecs$n.vars ){
		stop( 'Argument "sex" cannot be larger than "n.vars"', call. = F )
	}
}
#
if( F ){
	## CHECK THAT comb.sex IS NOT SPECIFIED UNLESS xchrom = T
	if( !.xchrom & !is.null( .info$model$comb.sex ) ) warning( 'Argument "comb.sex" is only implemented for models where "xchrom = TRUE"', call. = F )
}
#
## SET VARIABLE sel.sex SO THAT DON'T HAVE TO REPLACE WITH comb.sex ALL THE WAY THROUGH
if( identical( .info$model$comb.sex, "males" ) ) .info$variables$sel.sex <- 1
if( identical( .info$model$comb.sex, "females" ) ) .info$variables$sel.sex <- 2
#
## MAKE SURE FULL DATA'S NOT DUMPED WHEN scoretest = "only"
if( .info$model$scoretest == "only" & ( .info$control$data.out == "full" ) ){
	warning( 'Since data.out = "full", scoretest argument is changed from "only" to "no"', call. = F )
	.info$control$scoretest <- "no" 
}
#
## CHECK VALUES FOR MARKERS
if( !is.numeric( .info$filespecs$markers ) & !identical( .info$filespecs$markers, "ALL" ) ) stop( '"markers" argument must be either "ALL" ( default ) or an integer value.', call. = F ) 
#
## NO MATERNAL EFFECTS WITH PURE cc
if( ( .info$model$design == "cc" ) & ( .info$model$maternal ) ) stop( 'Cannot use maternal = TRUE with design = "cc"', call. = F )
#
## CHECK THAT response IS RESTRICTED
if( ( .info$haplos$response %in% c( "mult" ) ) & is.element( .info$haplos$reference, c( "reciprocal", "population" ) ) ){
	warning( 'response = "mult" must be used with reference category ( numeric or "ref.cat" ). Has been changed to reference = "ref.cat"', call. = F )
	.info$haplos$reference <- "ref.cat"
}
#
## IF test.maternal IS TRUE, MAKE SURE maternal IS ALSO TRUE
if( .info$model$test.maternal ) .info$model$maternal <- TRUE

## check name of given strata - should pick the env.exposure of the child!
if( !is.null( .info$variables$strata ) ){
	if( length( .info$variables$strata ) != 1 ){
		stop( "You gave more than one 'stratum'!", call. = FALSE )
	}
	if( !is.numeric( .info$variables$strata ) ){
		which.col.stratum <- match( paste0( .info$variables$strata, ".c" ), colnames( mcall$data$cov.data ) )
		if( is.na( which.col.stratum ) ){
			stop( "The given stratum name: ", .info$variables$strata, " doesn't match any of the covariate column names!", call. = FALSE )
		}
		.info$variables$strata <- which.col.stratum
	}
}

## If any of the parameters in 'info' were not set, there might be 'NAs' in the list, so let's remove it before returning
check.nas <- lapply( .info[ -(1:2) ], function(x){ any(is.na( names( x ) )) } )
if( any( unlist( check.nas ) ) ){
	cur.info <- .info[ -(1:2) ]
	which.nas <- which( unlist( check.nas ) )
	for( i in which.nas ){
		cur.list <- cur.info[[ i ]]
		which.nas2 <- which( is.na( names( cur.list ) ) )
		fixed.list <- cur.list[ -which.nas2 ]
		cur.info[[ i ]] <- fixed.list
	}
	.info[ -(1:2) ] <- cur.info
}

return( .info )


}
