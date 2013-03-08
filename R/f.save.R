f.save <- function(data, database, cpus = 1, batchsize = 1000, filehash = T, verbose = T){
##
## DUMP A FILE TO A DIRECTORY NAMED database, SPLIT OVER SEVERAL BATCHES TO AVOID
## LIMITATIONS ON NUMBER OF FILES AND MEMORY USAGE
##
#
### PREPARE
.para <- "snowfall"
#.para <- "doSMP" ## BURDE IMPLEMENTERES, SOM I haplinSlide?
#
.n <- dim(data)[2]
#
## GET BATCH SPLITUP
.batch <- f.batch(n = .n, batchsize = batchsize)
.ubatch <- unique(.batch$batch)
#
## FUNCTION TO DUMP EACH BATCH
.f.dump <- function(x){
	## DUMP INDATA FILE TO DATABASE DIRECTORY
	.cols <- .batch$col[.batch$batch == x]
	if(!filehash){
		.dir <- paste(database, "/batch_", x, ".RData", sep = "")
		if(verbose) cat("batch", x, "/", length(.ubatch), "\n")
		batch <- data[, .cols]
		save(batch, file = .dir)
	}else{
		library(filehash)## FOR THE SAKE OF snowfall AND doSMP
		.dir <- paste(database, "/batch_", x, sep = "")		
		if(verbose) cat("batch", x, "/", length(.ubatch), "\n")
		batch <- data[, .cols]
		dumpDF(batch, dbName = .dir, type = "RDS")
	}
}
#
## BEFORE DUMPING, CHECK NON-EXISTENCE
if(file.exists(database)) stop(paste("A directory, file, or database with the name ", database, " already exists!", sep = ""), call. = F)
dir.create(path = database)
#
## DUMP FILES
## EITHER IN SEQUENCE, SINGLE CPU
## OR IN PARALLEL (MULTIPLE CPUs, FASTER)
if(cpus == 1){
	lapply(.ubatch, .f.dump)
}else{
#	library(snowfall)
	if(!is.numeric(cpus)) stop('The number of cpu-s "cpus" must be numeric!', call. = F)
	#
	## INITIALIZE CPUS
	sfInit(parallel = T, cpus = cpus)
	on.exit(sfStop())
	#
	## RUN BY SPLITTING ON CPUS
	sfLapply(.ubatch, .f.dump)
}
#
## RETURN RELEVANT BATCH INFORMATION
return(.batch)
}

