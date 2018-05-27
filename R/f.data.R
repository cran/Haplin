f.data <- function( data.read, info, stratum.col, sel.stratum, quick = FALSE ){
##
## MISCELLANEOUS DATA PROCESSING
## 
##
## SET PARAMETERS, FOR SIMPLICITY
design <- info$model$design
xchrom <- info$model$xchrom
use.missing <- info$model$use.missing
verbose <- info$control$verbose
n.vars <-  info$filespecs$n.vars
markers <- info$filespecs$markers
family <- info$model$fam
sel.markers <- info$control$sel.markers

# choose the sex columns
.ind.sub <- TRUE
if( !missing( stratum.col ) ){
	.ind.sub <- ( data.read$cov.data[ ,stratum.col ] == sel.stratum )
}
if( xchrom & !is.null( info$variables$sel.sex ) ){
	## IF ON X-CHROM, AND ONLY ONE SEX IS SELECTED
	.sex <- data.read$cov.data[ ,info$variables$sex ]
	if( any( is.na( .sex ) ) ) {
		stop(paste0(sum(is.na(.sex)), " missing values found in sex variable! Must be removed from file before analysis.\n"), call. = F)
	}
	.tmp <- sort( unique( .sex ) )
	if( any( !is.element( .tmp, c("1", "2") ) ) ) {
		stop(paste0("The sex variable is coded ", paste(.tmp, collapse = " "), ". It should be coded 1 (males) and 2 (females). Missing values are not allowed."), call. = F)
	}

	if( length( .ind.sub ) != 1 ){
		.ind.sub <- ( .ind.sub | (.sex == info$variables$sel.sex) )
	} else {
	.ind.sub <- (.sex == info$variables$sel.sex)
	}
	## CHECK AND REPORT IF NO LINES LEFT
	if( sum( .ind.sub ) == 0 ){
		stop(paste0('No lines of data left when using "comb.sex = ', info$model$comb.sex, '"'), call. = F)
	}
}

## selecting the markers:
marker.names.all <- data.read$aux$marker.names
.markers <- TRUE
all.markers <- sum( sapply( data.read$gen.data, ncol ) )
if( sel.markers ){
	.sel <- f.sel.markers( n.vars = 0, markers = markers, family = family, split = T, ncols = all.markers )
	.markers <- attr( .sel, "markers" )
	info$filespecs$markers <- .markers
# 	info <- c( info, list( data = list( marker.names.all = marker.names.all[ .markers ] ) ) )
# 	info$data <- c( info$data, marker.names.all = marker.names.all[ .markers ] )
	info$data$marker.names.all <- marker.names.all[ .markers ]
} else {
	.sel <- 1:all.markers
}
#
## EXTRACT DATA COLUMNS AND ROWS
gen.data.extracted <- f.get.gen.data.cols( data.read$gen.data, as.numeric( .sel ) )
if( !all( .ind.sub ) ){
	gen.data.extracted <- f.extract.ff.numeric( gen.data.extracted, rows = which( .ind.sub ) )
} else {
	gen.data.extracted <- f.extract.ff.numeric( gen.data.extracted )
}
cov.data.extracted <- data.read$cov.data[ .ind.sub, ]
#
## HANDLE MISSING DATA:
.sum.na <- rowSums( is.na( gen.data.extracted ) )
.is.na <- .sum.na > 0.1
.rows.with.na <- sum( .is.na )
#
## REMOVE ROWS WITH MISSING, IF REQUESTED:
if( !use.missing ){
	if( all( .is.na ) ){
		stop('All data lines contain at least one missing, and "use.missing" is set to FALSE', call. = F)
	}
	.rows.dropped <- which( .is.na )
	
	if( n.vars > 0 ){
		data.new <- list( cov.data = cov.data.extracted[ !.is.na, ], gen.data = gen.data.extracted[ !.is.na, ] )
	} else {
		data.new <- list( cov.data = NULL, gen.data = gen.data.extracted[ !.is.na, ] )
	}
}else{
	.rows.dropped <- NULL
	
	data.new <- list( cov.data = cov.data.extracted, gen.data = gen.data.extracted )
}
# just to make sure that cov.data is a matrix
if( n.vars > 0 & is.null( dim( data.new$cov.data ) ) ){
	data.new$cov.data <- matrix( data.new$cov.data, ncol = 1 )
}

## check how many rows of data left
.nrow <- nrow( data.new$gen.data )
if( .nrow == 0 ){
	stop( "No rows of data available!", call. = F )
}
if( .nrow <= 15 ){
	warning( paste("Only ", .nrow, " lines of data available!\nHaplin is unlikely to produce very reliable results...." ), call. = F )
}


## COUNT AND REPORT MISSING
info$data$rows.dropped <- .rows.dropped

.marker.names <- marker.names.all[ .markers ]
data.new$aux$marker.names <- .marker.names

# THE NUMBER OF TRIADS AVAILABLE AT EACH STAGE:
.ntri.seq <- rep(NA, 4)
# THE ORIGINAL LINE NUMBERS AVAILABLE AT EACH STAGE:
.orig.lines.seq <- vector(4, mode = "list")
# NOTE: .ntri.seq[i] SHOULD BE THE SAME AS length(.orig.lines.seq[[i]]) FOR i = 1, 2
# WARNING: ALSO, .orig.lines.seq[[2]] SHOULD BE IN THE SAME ORDER AS THE CORRESPONDING DATA SET data.read!

names( .ntri.seq ) <- names( .orig.lines.seq ) <- c( "Original", "After rem NA", "After rem Mend. inc.", "After rem unused haplos" )

.ntri.seq[2] <- nrow( data.new$gen.data )

if(.rows.with.na == 0){
	.ntri.seq[1] <- .ntri.seq[2]
	.orig.lines.seq[[1]] <- .orig.lines.seq[[2]] <- 1:(.ntri.seq[1])
	if(verbose) cat("No lines contained missing data\n")
} else {
	if(use.missing){
		.ntri.seq[1] <- .ntri.seq[2]
		.orig.lines.seq[[1]] <- .orig.lines.seq[[2]] <- 1:(.ntri.seq[1])
		if(verbose) cat("There were ", .rows.with.na, " rows with missing data\nAll rows retained in analysis\n", sep = "")
	} else{
		.ntri.seq[1] <- .ntri.seq[2] + .rows.with.na
		.orig.lines.seq[[1]] <- .orig.lines.seq[[2]] <- 1:(.ntri.seq[1])
		.orig.lines.seq[[2]] <- .orig.lines.seq[[2]][-.rows.dropped]
		if(verbose){
			cat("The following", .rows.with.na, "data lines were dropped due to missing data:\n", .rows.dropped, "\n")
		}
	}
}

#------ done now in pre-processing -------
# ## FREQUENCY COUNT AND ALLELE SORTING:
# if(verbose) cat("\nPreparing data for analysis...  ")
# .data <- f.prep.data( data.read, info = info )
# if(verbose) cat("Done\n")

## EXTRACT ALLELE INFORMATION:
info$haplos$alleles <- data.read$aux$alleles[ .markers ]
info$haplos$alleles <- lapply( info$haplos$alleles, function(x){
	if( is.na( names(x)[ length(x) ] ) ){
		x[ -( length(x) ) ]
	} else {
		x
	}
} )
# NUMBER OF MISSING AT EACH LOCUS:
info$haplos$alleles.nas <- unlist( data.read$aux$alleles.nas[ .markers ] )

## CHANGE CASE: UPPER-CASE IS MOST FREQUENT
.f.change.case <- function(allele){
	names(allele) <- casefold(names(allele), upper = F)
	.max <- which(allele == max(allele))[1]
	names(allele)[.max] <- casefold(names(allele)[.max], upper = T)
	allele
}
info$haplos$alleles <- lapply( info$haplos$alleles, .f.change.case )

if( !is.null( .rows.dropped ) | !all( .ind.sub ) ){
	tmp.freq.list <- f.extract.freq( data.new, info$haplos$alleles, design )
	info$haplos$alleles <- tmp.freq.list$alleles
	info$haplos$alleles.nas <- unlist( tmp.freq.list$alleles.nas )
	data.new$aux$variables <- tmp.freq.list$variables
	data.new$aux$variables.nas <- tmp.freq.list$variables.nas
} else {
	data.new$aux$variables <- data.read$aux$variables
	data.new$aux$variables.nas <- data.read$aux$variables.nas
}

data.new$aux$alleles <- info$haplos$alleles
data.new$aux$alleles.nas <- info$haplos$alleles.nas

if( xchrom & !is.null( data.new$cov.data ) ){
	# need to re-code as numeric data
	data.new$cov.data[, info$variables$sex] <- as.numeric( data.new$cov.data[, info$variables$sex] )
}

## RETURN DATA BEFORE "HEAVY" PREPARATION
if(quick){
	data.out.tmp <- make.ff.data.out( covd = data.new$cov.data, gend = data.new$gen.data, covd.names = colnames( data.new$cov.data ), gend.names = colnames( data.new$gen.data ) )
	data.out.tmp$aux <- data.new$aux
	attributes( data.out.tmp ) <- attributes( data.new ) # xxx sjekk denne

	return(list( data = data.out.tmp, info = info ))
}

## DESIGN-DEPENDENT DATA PREPARATIONS:

## ORGANIZE GENETIC DATA,
## REMOVE MEND. INCONS.,
## ADD FREQUENCY COUNTER,
## SEPARATE INTO VARIABLES AND GENETIC DATA,
## (AND TEST FOR HWE)
# data.new$aux$variables <- data.new$aux$variables

.tmp <- f.sep.data( data.new, info )
.data.gen <- .tmp$data.gen
.HWE.res <- .tmp$HWE.res
.orig.lines.after.NA <- attr( .data.gen, "orig.lines" )
# A LIST OF THE ORIGINAL LINE NUMBERS (REFERS TO THE FILE AFTER POSSIBLE REMOVAL OF MISSING, THEN LINES WITH MEND.CONS. HAVE BEEN DELETED). CAN BE INDEXED BY ind.unique.line

## CHECK SOME OF THE HWE RESULTS:
if( !xchrom ){
	for( i in seq( along = info$haplos$alleles ) ){
		if( any( info$haplos$alleles[[i]] != .HWE.res[[i]]$freq ) ){
			warning( "Something's strange with the frequency count in HWE test!", call. = FALSE )
		}
	}
}

## REPORT MENDELIAN INCONSISTENCIES:
# LINE NUMBERS REFER TO DATA AFTER POSSIBLE REMOVAL OF MISSING 
.rows.with.Mendelian.inconsistency <- attr( .data.gen, "rows.with.Mendelian.inconsistency" )

if( length( .rows.with.Mendelian.inconsistency ) == 0 ){
	.ind.Mend <- numeric()
	.ntri.seq[3] <- .ntri.seq[2]
	.orig.lines.seq[[3]] <- .orig.lines.seq[[2]]
	if(!use.missing & .rows.with.na > 0){
		if(verbose) cat("None of the retained lines contained Mendelian inconsistencies\n")
	} else {
		if(verbose) cat("No lines contained Mendelian inconsistencies\n")
	}
}else{
	.ind.Mend <- .orig.lines.seq[[2]][.rows.with.Mendelian.inconsistency] ## WILL REFER TO LINE NUMBERS (WITH POSS. MEND. INCONS.) IN ORIGINAL FILE
	if(verbose){
		cat( "The following", length(.ind.Mend), "data lines were dropped due to Mendelian inconsistencies:\n", .ind.Mend, "\n" )
	}
	.orig.lines.seq[[3]] <- .orig.lines.seq[[2]][ -.rows.with.Mendelian.inconsistency ]
	.ntri.seq[3] <- length( .orig.lines.seq[[3]] )
}


## PRELIMINARY DATA FIXUP:
## EXPAND FREQUENCIES AND ADD COUNTER:
.orig.lines.after.NA <- unlist( .orig.lines.after.NA[ .data.gen$ind.unique.line ] )
# CONVERT LINE NUMBERS INTO THE ORIGINAL LINE NUMBERS
.orig.lines <- .orig.lines.seq[[2]][.orig.lines.after.NA]
# WARNING: .orig.lines.seq[[2]] SHOULD HAVE THE SAME ORDERING AS data.read!
if( any( .orig.lines.seq[[3]] != sort( unique( .orig.lines ) ) ) ){
	stop( "problem!", call. = FALSE )
}

.ind <- 1:( dim(.data.gen)[1] )
.ind <- rep( .ind, .data.gen$freq )
# .ind.aux <- unlist( sapply( .data.gen$freq, function(x) 1:x ) )


if(design == "triad" | design == "cc.triad"){
	if(!xchrom){
		col.idx <- 1:5
		new.data.gen.names <- c( "m1", "m2", "f1", "f2", "ind.unique.line", "orig.lines" )
	}
	if(xchrom){
		col.idx <- 1:6
		new.data.gen.names <- c( "m1", "m2", "f1", "f2", "sex", "ind.unique.line", "orig.lines" )
	}
}
if(design == "cc"){
	if(!xchrom){
		col.idx <- 1:3
		new.data.gen.names <- c( "c1", "c2", "ind.unique.line", "orig.lines" )
	}
	if(xchrom){
		col.idx <- 1:4
		new.data.gen.names <- c( "c1", "c2", "sex", "ind.unique.line", "orig.lines" )
	}
}

ncol.idx <- length( col.idx )
.data <- dframe( .data.gen[ .ind,col.idx ], .orig.lines )
names( .data ) <- new.data.gen.names

## REPLACE LINE COUNTERS ETC. WITH UNIQUE TAG ind, WHICH HAS ONE VALUE FOR EACH (REMAINING) TRIAD:
.data$ind.unique.line <- NULL

## COMPUTE PRELIMINARY HAPLOTYPE FREQUENCIES USING A SIMPLE EM-VERSION:
.prelim.freq <- f.preliminary.freq.new( .data, info )
.data$freq <- ff::as.ff( .prelim.freq )
info$haplos$prelim.haplotype.freq <- attr( .prelim.freq, "prelim.haplotype.freq" )

## DECIDE WHICH HAPLOTYPES TO INCLUDE IN ANALYSIS
info$haplos$selected.haplotypes <- f.sel.haplos( info )
.n.sel.haplos <- sum( info$haplos$selected.haplotypes )

## REMOVE HAPLOTYPES WITH INITIAL FREQUENCY BELOW threshold.
## HAPLOTYPES ARE RECODED TO 1,2,3,... AFTER REMOVAL.
## FREQUENCIES ARE RENORMALIZED SO THAT EACH TRIAD SUMS TO ONE.
##
if(verbose) cat( "Removing unused haplotypes...  " )

if(abs(sum(.data$freq[]) - .ntri.seq[3]) > 1e-6){
	warning( "There may be a problem with the data summary", call. = FALSE )
}
.data <- f.repl.thin(.data, selection = info$haplos$selected.haplotypes, design = design)
attr( .data, "selected.haplotypes" ) <- info$haplos$selected.haplotypes
.ntri.seq[4] <- sum(.data$freq[])
if(abs(.ntri.seq[4] - round(.ntri.seq[4])) > 1e-6){
	warning( "There may be a problem with the data summary", call. = FALSE )
}
.ntri.seq[4] <- round(.ntri.seq[4])
.orig.lines.seq[[4]] <- unique(.data$orig.lines[])
if(verbose) cat("Done\n")

## DECIDE REFERENCE
.tmp <- f.prep.reference(info)
info$haplos$reference.method <- .tmp$reference.method
info$haplos$ref.cat <- .tmp$ref.cat

## ADD ON CASE-CONTROL VARIABLE FOR cc AND cc.triad DATA
if(design == "cc" | design == "cc.triad"){
	.ccvar <- info$variables$ccvar
	.tmpind <- match(.data$orig.lines[], .orig.lines.seq[[2]])
## WARNING: .data.vars SHOULD STILL HAVE THE SAME ORDERING AS data.read, AND .orig.lines.seq[[2]] SHOULD REFER TO THIS ORDERING!
	.cc <- as.numeric( unlist( data.new$cov.data[.tmpind, .ccvar] ) )
	if(any(is.na(.cc))) {
		stop(paste(sum(is.na(.cc)), " missing values found in case-control variable! Must be removed from file before analysis.\n", sep = ""), call. = F)
	}
	.cctab <- data.new$aux$variables[[ .ccvar ]]
	.codes <- names( .cctab )
	if(length(.codes) != 2) {
		stop(paste('Case-control variable "ccvar" is coded as ', paste(.codes, collapse = ", "), '. It should have exactly two different values!', sep = ""), call. = F)
	}
	if(!identical(sort(unique(.cc)), c(1,2))) {
		stop("Something's wrong with the case-control variable!", call. = F)
	}
	if(verbose){
		cat("\nNote:\n")
		cat("The selected case/control variable is column ", .ccvar, ": '", colnames(data.new$cov.data)[.ccvar], "'\n", sep = "")
		cat("The following case/control coding and frequencies have been found:\n")
		cat("controls(", .codes[1], "): ", .cctab[1], ", cases(", .codes[2], "): ", .cctab[2], "\n", sep = "")
	}
	.data$cc <- .cc
}

## ADD ON COVARIATE INFORMATION IF REQUESTED
if(!is.null(info$variables$covar)){
	stop('The "covar" argument is not available in "haplin" and "haplinSlide", only in "haplinStrat"!', call. = F)
	.covar <- info$variables$covar ## COLUMN NUMBER
	.tmpind <- match(.data$orig.lines[], .orig.lines.seq[[2]])
	## WARNING: .data.vars SHOULD STILL HAVE THE SAME ORDERING AS data.read, AND .orig.lines.seq[[2]] SHOULD REFER TO THIS ORDERING!
	.co <- as.numeric( unlist( data.read$cov.data[.tmpind, .covar] ) ) ## INTEGER VALUES FROM RECODED DATA FILE
	if(any(is.na(.co))) {
		stop(paste(sum(is.na(.co)), " missing values found in covariate! Must be removed from file before analysis.\n", sep = ""), call. = F)
	}
	## COVAR TABLE
	.covar.tab <- data.new$aux$variables[[.covar]]
	if(length(.covar.tab) == 1) {
		stop(paste('Covariate variable "covar" has only a single value ', paste(.codes, collapse = ", "), '. It should have two or more different values!', sep = ""), call. = F)
	}
	## (A PROBABLY UNNECESSARY) QUICK CHECK:
	.tmpco <- sort(unique(.co))
	if(any(.tmpco != seq(along = .tmpco))) {
		stop('Something\'s wrong with the "covar" variable!', call. = F)
	}
	
	if(verbose){
		cat("\nDistribution of the covariate variable:", sep = "")
		print(data.new$aux$variables[[.covar]])
	}
	.data$covar <- .co
	## ADD COVARIATE INFORMATION TO .info
	info$variables$covar.codes <- .covar.tab
}

## ADD INFORMATION TO THE .info OBJECT
info$data$ntri.seq <- .ntri.seq
info$data$lines.Mend.inc <- .ind.Mend ## LINE NUMBERS (IN ORIGINAL FILE) WITH MEND. INCONS.
info$check$HWE.res <- .HWE.res

return( list( data = .data, info = info ) )
}
