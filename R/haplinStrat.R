haplinStrat <- function( data, strata = NULL, ... ){
##
## RUN HAPLIN _FIRST_ ON FULL DATA SET, 
## THEN, INDEPENDENTLY, ON EACH STRATUM DEFINED BY strata.
#
## GET HAPLIN DEFAULTS, OVERRIDE DEFAULTS FOR SPECIFIED ARGUMENTS
.info <- f.catch(match.call(), formals())
.verbose <- .info$control$verbose
.use.missing <- .info$model$use.missing

if(.verbose){
	cat("\n## Running haplinStrat ##\n")
}

if( !is.null( attr( data, "sel.markers" ) ) ){
	.info$control$sel.markers <- attr( data, "sel.markers" )
}

if( !(is(data, "haplin.ready") | is(data, "prep.data")) ){
	stop( "The given input data is not ready for haplin analysis! Please use 'genDataPreprocess' function first!", call. = FALSE )
}

all.markers <- sum( unlist( lapply( data$gen.data, ncol ) ) )

### DEFINE STRATA VARIABLE AND DISPLAY FREQUENCY DISTRIBUTION
.strata <- data$cov.data[, .info$variables$strata]
.strata.names <- names(data$aux$variables[[.info$variables$strata]])
#
if( !.use.missing & (data$aux$nrows.with.missing > 0) ){
	.strata <- .strata[ -data$aux$which.rows.with.missing ]
} 
#
if(anyNA(.strata)) stop("Missing values found in stratification variable.\n\ \ Should be removed from file before running haplinStrat.", call. = F)

if(.verbose){
	cat("\nSelected stratification variable:", names(data$cov.data)[.info$variables$strata])
	cat("\nFrequency distribution of stratification variable:\n")
	.tmp.tab <- table(.strata.names[.strata])
	names(dimnames(.tmp.tab)) <- NULL
	print(.tmp.tab)
}
#
## PREPARE RESULTS LIST
.strata.list <- sort(unique(.strata))
.ut.list <- vector(length(.strata.list) + 1, mode = "list")
names(.ut.list) <- c("all", .strata.names[.strata.list])
#
## SET UP TEMPORARY FILE FOR HAPLOTYPES
.tmphaplofile <- tempfile(tmpdir = ".")
on.exit(unlink(.tmphaplofile))
#
## PREPARE ARGUMENTS
.args <- f.args.from.info(.info)
# .args$filename <- NULL # REMOVE ANY FILENAME, FROM NOW ON USE ONLY .data
all.column.names <- lapply(
  data$gen.data,
  colnames
) |>
  do.call(
    what = c, args = _
  )

if( .info$control$sel.markers ){
	.sel <- f.sel.markers(
	  n.vars = 0,
	  markers = .info$filespecs$markers,
	  family = .info$model$fam,
	  split = T,
	  all.marker.names = all.column.names
	)
	cur.markers <- attr( .sel, "markers" )
	gen.data.tmp <- f.get.gen.data.cols( data$gen.data, .sel )
	if( !is.null( data$cov.data ) ){
		data.tmp <- make.ff.data.out( covd = data$cov.data, gend = list( gen.data.tmp ), data.as.is = TRUE )
	} else {
		data.tmp <- make.ff.data.out( gend = list( gen.data.tmp ), data.as.is = TRUE )
	}
	data.tmp$aux <- data$aux # neppe god
	attributes( data.tmp ) <- attributes( data )# xxx trengs?
	# attr( data.tmp, "alleles" ) <- attr( data, "alleles" )[ cur.markers ]
	# attr( data.tmp, "alleles.nas" ) <- attr( data, "alleles.nas" )[ cur.markers ]
	# attr( data.tmp, "marker.names" ) <- attr( data, "marker.names" )[ cur.markers ]

	data.tmp$aux$alleles <- data$aux$alleles[ cur.markers ]
	data.tmp$aux$alleles.nas <- data$aux$alleles.nas[ cur.markers ]
	data.tmp$aux$marker.names <- data$aux$marker.names[ cur.markers ]
} else {
	data.tmp <- data
}
.info$control$sel.markers <- FALSE

.args$markers <- "ALL" # HAVE ALREADY SELECTED RELEVANT MARKERS
attr( data.tmp, "sel.markers" ) <- FALSE

.args$strata <- NULL # SHOULD NOT BE SENT TO haplin
.args$verbose <- F # SHOULD PERHAPS ALLOW MORE FLEXIBILITY HERE?
.args$printout <- F # SHOULD PERHAPS ALLOW MORE FLEXIBILITY HERE?
#
## RUN ON FULL DATA
if(.verbose) cat("\nRunning Haplin on full data file...")
#
## SET DATA TO FULL FILE
.args$data <- data.tmp

## AND RUN
.ut.list[["all"]] <- do.call("haplin", .args)
if(.verbose) cat("Done\n")
#
## WRITE (TEMPORARY) FILE CONTAINING HAPLOTYPES
.selected.haplotypes <- .ut.list[["all"]]$info$haplos$selected.haplotypes
.selected.haplotypes <- names(.selected.haplotypes)[.selected.haplotypes]
write.table(dframe(haplos = .selected.haplotypes), file = .tmphaplofile, quote = F, row.names = F, col.names = T)
## FORCE haplin LATER ON TO USE SAME HAPLOTYPES AND SAME REFERENCE CATEGORY
.args$haplo.file <- .tmphaplofile
.args$reference <- .ut.list[["all"]]$info$haplos$ref.cat
.reference.method <- .ut.list[["all"]]$info$haplos$reference.method
#
## RUN HAPLIN ON EACH STRATUM
for(i in seq(along = .strata.list)){
	.mess <- paste('\nRunning Haplin on stratum "', .strata.names[.strata.list[i]], '"...', sep = "") ## NEED TO DEF. THIS MESSAGE BEFOREHAND, OTHERWISE cat DOESN'T PRINT EVERYTHING AT ONCE!
	if(.verbose) cat(.mess)
	#
	## FEED haplin WITH STRATA SUBSET OF FILE
	tmp.info <- .ut.list[["all"]]$info
	tmp.info$control$sel.markers <- FALSE
	data.info.subset <- f.data( data.tmp, tmp.info , stratum.col = .info$variables$strata, sel.stratum = .strata.list[i], quick = TRUE )

	data.subset <- data.info.subset$data
	
	data.subset$gen.data <- list( data.subset$gen.data ) # xxx why??
	class( data.subset ) <- "haplin.ready"
	attr( data.subset, "sel.markers" ) <- FALSE
	.args$data <- data.subset

	## RUN HAPLIN
	.ut.list[[i + 1]] <- try(do.call("haplin", .args))
	if(!is(.ut.list[[i + 1]], "try-error")){
		## RESET reference.method TO ORIGINAL
		if(.ut.list[[i + 1]]$info$haplos$reference.method != "ref.cat") stop("Something's wrong with the choice of reference...", call. = F) # should never kick in
		.ut.list[[i + 1]]$info$haplos$reference.method <- .reference.method
	}
	#
	if(.verbose) cat("Done\n")
}
#
class(.ut.list) <- "haplinStrat"
#
return(.ut.list)
}
