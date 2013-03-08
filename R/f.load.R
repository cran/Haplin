f.load <- function(database, n, batchsize = 1000, cols = 1:3, varnames = c("x1", "x2", "x3"), filehash = T){
##
## SELECT AND LOAD COLUMNS FROM A "HAPLIN" DATABASE
##
## database IS THE (CHARACTER) NAME OF DIRECTORY CONTAINING THE DATABASE
## n IS THE TOTAL NUMBER OF COLUMNS IN DATABASE (SHOULD BE IMPROVED!)
## batchsize IS THE BATCH SIZE THE DATABASE WAS SAVED WITH
## cols ARE THE COLUMNS TO BE PICKED
## varnames ARE THE CORRESPONDING VARIABLE NAMES
##
##

# n KUNNE VAERE STORESTE VERDI AV cols?

## FIND WHAT BATCHES THE COLUMNS BELONG TO
.batch <- f.batch(n = n, batchsize = batchsize)
#
## SELECT ONLY RELEVANT BATCHES/COLUMNS
.batch <- .batch[.batch$col %in% cols,]## SJEKK HVA SOM SKJER HVIS GAAR UTENFOR
#
## UNIQUE BATCHES
.ubatch <- unique(.batch$batch)
#
## PREPARE TO SELECT
.res.list <- vector(length(.ubatch), mode = "list")
for(i in seq(along = .res.list)){
	if(!filehash){
		## BATCH FILE NAME
		.dir <- paste(database, "/batch_", .ubatch[i], ".RData", sep = "")
		## RELEVANT COLUMN NAMES FOR THIS BATCH
		.cols <- .batch$col[.batch$batch == .ubatch[i]]
		.navn <- varnames[.cols]
		###.navn <- paste("V", .cols, sep = "")
		###.db <- dbInit(.dir, type = "RDS")
		###.res.list[[i]] <- dbMultiFetch(.db, .navn)
		## LOAD DATA LOCALLY, EXTRACT COLUMNS
		###.res.list[[i]] <- local({ ### HAR ERSTATTET MED DET NEDENFOR
		###	load(.dir)
		###	batch[, .navn]
		###})
		.res.list[[i]] <- within(list(), load(.dir))$batch[,.navn] ### HAR IKKE SJEKKET....
	}else{
		## BATCH FILE NAME
		.dir <- paste(database, "/batch_", .ubatch[i], sep = "")
		## RELEVANT COLUMN NAMES FOR THIS BATCH
		.cols <- .batch$col[.batch$batch == .ubatch[i]]
		.navn <- varnames[.cols]
		.tmp <- dbInit(.dir, type = "RDS")
		.res.list[[i]] <- dbMultiFetch(.tmp, .navn)
	}
}
#
## CONVERT TO DATA FRAME
.ut <- as.dframe(.res.list)
#
##
return(.ut)
}


