haplinStrat <- function(filename, data, pedIndex, covar = NULL, ...){
#
#
## GET HAPLIN DEFAULTS, OVERRIDE DEFAULTS FOR SPECIFIED ARGUMENTS
.info <- f.catch(match.call(), formals())


.verbose <- .info$control$verbose
#
##
if(.verbose){
	cat("\n## Running haplinStrat ##\n")
	cat("\nReading data from file...  ")
}
#
## GET FULL DATA
.data.read <- f.get.data(data, pedIndex, .info)
## UPDATE .info
.info <- attr(.data.read, "info")


#
### DEFINE STRAT VARIABLE AND DISPLAY FREQUENCY DISTRIBUTION
.strata <- .data.read[, .info$variables$covar]
if(.verbose){
	cat("\nFrequency distribution of selected stratification variable:\n")
	.tmp.tab <- table(.strata)
	names(dimnames(.tmp.tab)) <- NULL
	print(.tmp.tab)
}
if(any(.strata == "NA")) stop("Missing values found in stratification variable.\n\ \ Should be removed from file before running haplinStrat.")
### MERK, MERK: I f.read.data BLIR MISSING KONVERTERT TIL NA OG DERETTER "NA". BURDE KANSKJE BARE VAERT NA?
#
## PREPARE RESULTS LIST
.strata.list <- sort(unique(.strata))
.ut.list <- vector(length(.strata.list) + 1, mode = "list")
names(.ut.list) <- c("all", .strata.list)
#
## SET UP TEMPORARY FILES
.tmpfilename <- tempfile(tmpdir = ".")
.tmphaplofile <- tempfile(tmpdir = ".")
on.exit({unlink(.tmpfilename); unlink(.tmphaplofile)})
#
## WRITE FULL DATA
write.table(.data.read, file = .tmpfilename, sep = "\t", col.names = F, row.names = F, quote = F)
#
## RUN ON FULL DATA
if(.verbose) cat("\nRunning Haplin for full data file...")

.args <- f.args.from.info(.info)


.argstmp <- .args
.argstmp$filename <- .tmpfilename
.argstmp$markers <- "ALL"
.argstmp$sep <- "\t"
.argstmp$allele.sep <- "\t"
.argstmp$na.strings <- "NA"
.argstmp$verbose <- F
.argstmp$printout <- F

.argstmp$covar <- NULL

.ut.list[["all"]] <- do.call("haplin", .argstmp)
if(.verbose) cat("Done\n")

#
## WRITE (TEMPORARY) FILE CONTAINING HAPLOTYPES
write.table(haptable(.ut.list[["all"]]), file = .tmphaplofile, quote = F)

.argstmp$haplo.file <- .tmphaplofile
.argstmp$reference <- .ut.list[["all"]]$info$haplos$ref.cat

#
## RUN HAPLIN FOR EACH STRATUM
for(i in seq(along = .strata.list)){
	.mess <- paste('\nRunning Haplin for stratum "', .strata.list[i], '"...', sep = "")
	if(.verbose) cat(.mess)
	.tmpd <- .data.read[.strata == .strata.list[i], ]
	## WRITE TEMPORARY HAPLIN DATA FILE
	write.table(.tmpd, file = .tmpfilename, sep = "\t", col.names = F, row.names = F, quote = F)
	## RUN HAPLIN
	.ut.list[[i + 1]] <- try(do.call("haplin", .argstmp))
	if(.verbose) cat("Done\n")
}

class(.ut.list) <- "haplinStrat"
return(.ut.list)


}
