toDataFrame <- function(x, reduce = F){
##
## STACK A LIST OF DATA FRAMES
.x <- x
#
## PREPARE MARKER NAMES
.markers <- names(.x)
.markers <- paste(.markers, "-xXx-", sep = "")
#
## FIRST NON-MISSING:
.is.na <- is.na(.x)
.start <- which(!.is.na)[1]
.first <- .x[[.start]]
.dim <- dim(.first)
.colnames <- colnames(.first)
#
## FILL IN THOSE THAT ARE MISSING, SO THAT THEY ARE RETAINED IN TABLE
## USE .first AS TEMPLATE
.miss <- .first[1,]
.miss[,] <- NA
.x[.is.na] <- replicate(sum(.is.na), .miss, simplify = F)
#
## UNLIST ONE LEVEL
.x <- unlist(.x, recursive = F)
#
## SET DIMENSION OF UNLISTED OBJECT
## NB: ASSUMES SAME NUMBER OF COLUMNS IN ALL ELEMENTS OF ORIGINAL LIST!
dim(.x) <- c(.dim[2], length(.x)/.dim[2])
colnames(.x) <- .markers
#
## STACK EACH COLUMN INDIVIDUALLY, USING UNLIST
.ut <- vector(.dim[2], mode = "list")
names(.ut) <- .colnames
for(i in 1:.dim[2]){
	.ut[[i]] <- unlist(.x[i,])
}
#
## SET IN DATA FRAME
.ut <- do.call("dframe", .ut)
#.ut <- as.dframe(.ut)
.markers.ext <- strsplit(rownames(.ut), split = "-xXx-")
.markers.ext1 <- sapply(.markers.ext, function(x)x[1])
.markers.ext2 <- sapply(.markers.ext, function(x)x[2])
.ut <- cbind(marker = .markers.ext1, row.no = .markers.ext2, .ut)
rownames(.ut) <- NULL
#
if(reduce){
	# REDUCE TO ONLY ONE LINE PR SNP IF ALL MARKERS ARE JUST SINGLE SNPS
	#
	## REMOVE LINE NUMBER IN THIS CASE
	.ut$row.no <- NULL
	## SELECT RELEVANT
	.row1 <- .ut[!duplicated(.ut$marker),]
	.ind.nonref <- is.na(.ut$reference) | .ut$reference != "ref"
	.row.nonref <- .ut[.ind.nonref,]
	if(any(.row.nonref$marker != .row1$marker)) stop()
	#c("marker", "Original", "After.rem.NA", "After.rem.Mend.inc.", "After.rem.unused.haplos", "alleles", "counts", "HWE.pv")
	.ut <- cbind(.row1[, 1:8], .row.nonref[, -(1:8)])
}
#
##
return(.ut)
}
