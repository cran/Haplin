f.batch <- function(n, batchsize, slice = 1){
##
## n IS THE NUMBER OF DATA COLUMNS
## batchsize IS THE NUMBER OF COLUMNS TO BE PUT IN EACH BATCH
## slice IS THE NUMBER OF COLUMNS BELONGING TO A SINGLE MARKER,
## MAKES SURE THEY'RE IN THE SAME BATCH (EVER REALLY USED?)
#
## 
.nslices <- n/slice
if(.nslices != round(.nslices)) stop('Something is wrong with the column number, or "slice" argument')
#
## NUMBER OF COMPLETE BATCHES
.hele <- .nslices %/% batchsize
#
## REMAINING COLUMNS
.rest <- .nslices - batchsize * .hele
#
##
.batch <- rep(seq(length.out = .hele), each = batchsize)
.batch <- c(.batch, rep(.hele + 1, length.out = .rest))
.batch <- rep(.batch, each = slice)
#
##
.slice <- rep(1:.nslices, each = slice)
#
##
if(F){
	#
	##
	.run <- rep(1:batchsize, each = slice)
	.run <- rep(.run, length.out = n)
	.reg <- dframe(col = seq(length.out = n), slice = .slice, batch = .batch, run = .run)
}else{
	.reg <- dframe(col = seq(length.out = n), slice = .slice, batch = .batch)
}
#
##
return(.reg)
}
