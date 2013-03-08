saveData <- function(filename, database, 
markers = "ALL", n.vars = 0, sep = " ", allele.sep = ";", na.strings = "NA", 
design = "triad", 
ccvar = NULL,
verbose = TRUE,
cpus = 1
){
##
## READS RAW DATA IN HAPLIN FORMAT. PREPARES AN info OBJECT FROM CALL
## USES f.read.data TO READ AND PROCESS. 
## THEN STORES DATA IN A "HAPLIN" DATABASE
##
## TO BE CALLED USING ARGUMENTS IDENTICAL TO HAPLIN'S FILE SPECIFICATIONS
##
## NOTE: use.missing IS ALWAYS SET TO TRUE HERE 
## USER CAN SELECT MARKERS, BUT NOT USE sel.sex OR SET use.missing  = F


# PROBLEM: HVIS MAN GLEMMER n.vars FAAR MAN BESKJED OM FEIL I allele.sep
# SJEKK HVIS DATABASE ER F.EKS. GUGG/


## CREATE info OBJECT
.defaults <- formals()
## AVOID WARNINGS FROM f.check.pars WHEN design == "cc"
.defaults$response <- "mult"
.defaults$reference <- "ref.cat"
## OVERRIDE STANDARD HAPLIN DEFAULT
.defaults$use.missing <- T
##
.info <- f.catch(match.call(), .defaults)
#.info <- f.catch(match.call(), formals())
#
## CHECK THAT filename IS NOT ALREADY A DATABASE
if(.info$filespecs$database){
	stop(paste(filename, " is already a database.", sep = ""))
}
#
## CHECK SOUNDNESS OF database ARGUMENT (WHICH IS NOT IN .info)
if(missing(database)) stop('File name argument "database" must be specified!', call. = F)
if(file.exists(database)) stop(paste("A directory, file, or database with the name ", database, " already exists!", sep = ""), call. = F)
#
## READ FULL DATA FILE
if(.info$control$verbose)	cat("\nReading data from file...  ")
.data.read <- f.read.data(info = .info) 
#
##
.info <- attr(.data.read, "info")
.data.read <- as.dframe(.data.read)
if(.info$control$verbose){
	cat("Done\n")
	cat("\nCreating data storage...\n")
}
#
## DUMP INDATA FILE TO DATABASE DIRECTORY
.dumpinfo <- f.save(data = .data.read, database = database, cpus = cpus, verbose = .info$control$verbose)
#
## SAVE RELEVANT DATABASE INFO:
#
## GIVE FILENAMES
.file.varnames <- paste(database, "/000varnames.RData", sep = "")
.file.info <- paste(database, "/000info.RData", sep = "")
.file.dumpinfo <- paste(database, "/000dumpinfo.RData", sep = "")
#
## SET REASONABLE OBJECT NAMES FOR SAVING
varnames <- names(.data.read)
info <- .info
dumpinfo <- .dumpinfo
#
## SAVE
save(varnames, file = .file.varnames)
save(info, file = .file.info)
save(dumpinfo, file = .file.dumpinfo)
if(.info$control$verbose) cat("Done\n")
#
##
return(invisible(.info))
}

