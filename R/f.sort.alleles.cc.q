f.sort.alleles.cc <- function(data)
{
## 
##
## PREPARATIONS:
	.alleles <- attr(data, "alleles")
	.nalleles <- sapply(.alleles, length)
	.nmarkers <- dim(data)[2]/2
#
#
## SET IN "CANONICAL" ORDERING, WITH SMALLEST FIRST:
#
	for(i in seq(1, .nmarkers * 2, 2)){
		.tmpcol1 <- pmin(data[,i], data[,i+1])
		.tmpcol2 <- pmax(data[,i], data[,i+1])
		data[,i] <- .tmpcol1	
		data[,i+1] <- .tmpcol2		
	}
#
## AGGREGATE DATA WITH A FREQUENCY COUNT:
	.data.agg <- f.aggregate(data)
	.lines <- attr(.data.agg, "orig.lines")
	.nlines.unique <- dim(.data.agg)[1]
#
#
## RESHAPE AND PERMUTE COLUMNS (ALL POSSIBLE PERMUTATIONS) WITHIN CHILD, RESULTING
## IN 2 COLUMNS (ROW)SORTED BY MARKER AND BY PERMUTATION WITHIN MARKER, 2 ROWS PR.
## CHILD PR MARKER.
#
	.data.long <- as.matrix(.data.agg[,1:(2 * .nmarkers)])
	.data.long <- t(matrix(as.numeric(t(.data.long)), nrow = 2)) # STACK DATA LINE BY LINE
	.swap <- c(c(1,2), c(2,1))
	.data.long <- .data.long[,.swap]
	.data.long <- t(matrix(as.numeric(t(.data.long)), nrow = 2))
	dimnames(.data.long) <- list(NULL, c("c1", "c2"))

## ADD A VECTOR ind.unique.line WHICH COUNTS LINE NUMBER AMONG THE UNIQUE LINES
## (DOES NOT REFER TO THE ORIGINAL LINE NUMBERS), AND A VECTOR ind.unique.line.marker
## WHICH IS A UNIQUE TAG FOR EACH UNIQUE LINE/MARKER COMBINATION:
	.data.long <- cbind(.data.long, ind.unique.line = rep(1:.nlines.unique, each = 2*.nmarkers), ind.marker = rep(1:.nmarkers, rep(2, .nmarkers))) 
###	.ind.unique.line.marker <- f.pos.in.grid(A = c(.nmarkers, .nlines), comb = .data.long[,c("ind.marker", "ind.unique.line")])
	.ind.unique.line.marker <- f.pos.in.grid(A = c(.nmarkers, .nlines.unique), comb = .data.long[,c("ind.marker", "ind.unique.line")])	
	.data.long <- cbind(.data.long, ind.unique.line.marker = .ind.unique.line.marker)
#
if(F){
#
## IDENTIFY MENDELIANLY CONSISTENT COMBINATIONS:
#
	.valid2 <- (.data.long[,2] == .data.long[,5]) # SECOND IN MOTHER EQUALS FIRST IN CHILD
	.valid2[is.na(.valid2)] <- T # IF ONE OR BOTH ARE MISSING THE COMBINATION IS VALID
	.valid4 <- (.data.long[,4] == .data.long[,6]) # SECOND IN FATHER EQUALS SECOND IN CHILD
	.valid4[is.na(.valid4)] <- T # IF ONE OR BOTH ARE MISSING THE COMBINATION IS VALID
#	
	.valid <- .valid2 & .valid4
	.valid.markers <- tapply(.valid, as.data.frame(.data.long[,c("ind.unique.line", "ind.marker")]), any)
	.valid.unique.lines <- apply(.valid.markers, 1, all) # NOTE: ASSUMES COMPLETE LIST OF INTEGERS FOR ind.unique.line, AND SORTED..
###	.rows.with.Mendelian.inconsistency <- unlist(.lines[.tag[.unique][!.valid.unique.lines]])	
	.rows.with.Mendelian.inconsistency <- unlist(.lines[!.valid.unique.lines])	
#
## REMOVE ALL INCONSISTENT LINES (WARNING: COMPLETELY INCONSISTENT WILL HERE DISAPP.!):
#
##	.data.long <- data.frame(.data.long, .valid)
	.keep <- .valid.unique.lines[.data.long[,"ind.unique.line"]] & .valid # KEEP ALL COSISTENT FOR THOSE WITH AT LEAST ONE CONSISTENT
	.data.long <- .data.long[.keep,]
#
## REDUCE TO A 4-COLUMN MATRIX (STILL LONG FORMAT) FOR ONLY MOTHER AND FATHER, REPLACE THE NAs WHEREEVER
# POSSIBLE, AND LEAVE REAL NAs OPEN:
#
	.data.long[,2] <- ifelse(!is.na(.data.long[,2]), .data.long[,2], .data.long[,5])
	.data.long[,4] <- ifelse(!is.na(.data.long[,4]), .data.long[,4], .data.long[,6])
#
	.data.long <- .data.long[,-c(5,6)]
}
#
## REMOVE DUPLICATE ROWS (FOR CC THESE ARE JUST ONE OF THE TWO HOMOZYGOTES):
#
	.tag.allcol <- f.create.tag(.data.long)
###	.unique.long <- unlist(tapply(.tag.4col, .data.long[,"ind.unique.line.marker"], function(x) !duplicated(x)))
	.ind.unique.allcol <- !duplicated(.tag.allcol)
	.data.long <- .data.long[.ind.unique.allcol,]
#
##	.l <- dim(.data.long)[1]
###	.nalleles.long <- .nalleles[.data.long[,"ind.marker"]]
#
## EXPAND ALL NAs WITH SEQUENCE OF ALL POSSIBLE ALLELES FOR THE CORRESPONDING MARKER:
#
	for(i in 2:1){# FOR EACH OF THE 2 ALLELES
		for(j in 1:.nmarkers){
			# HOW MANY LINES TO EXPAND EACH MISSING INTO:
			.nexpand <- rep(1, dim(.data.long)[1])			
			.nexpand[is.na(.data.long[,i]) & (.data.long[,"ind.marker"] == j)] <- .nalleles[j]
			# CREATE THE INDEX FOR EXPANSION AND EXPAND:
			.ind.expand <- rep(seq(along = .nexpand), .nexpand)
			.data.long <- .data.long[.ind.expand,]
			# REPLACE THE NAs IN THE EXPANDED DATA WITH SEQUENCE OF ALL POSSIBLE ALLELES AT MARKER:
			.data.long[is.na(.data.long[,i]) & (.data.long[,"ind.marker"] == j),i] <- 1:.nalleles[j]	
		}
	}
	# THIS MATRIX DOES NOT HAVE ANY REDUNDANT ROWS (?), BUT MAYBE AFTER PLACING ON SINGLE LINE...(?)
##
## SET MARKERS SIDE BY SIDE, WITH ALL POSSIBLE COMBINATIONS:
	.line.long <- 1:(dim(.data.long)[1])
	# MATRIX OF LINE NUMBERS CORRESPONDING TO EACH COMB OF UNIQUE LINES AND MARKERS:
	.line.bits <- tapply(.line.long, as.data.frame(.data.long[,c("ind.unique.line", "ind.marker")]), function(x)x)
	# CREATE ALL POSSIBLE COMBINATIONS OF MARKER GENOTYPES WITH OTHER MARKER GENOTYPES:
	.line.seq <- apply(.line.bits, 1, function(x) as.numeric(t(as.matrix(do.call("expand.grid", x)))))
	.line.seq <- unlist(.line.seq)	
	# EXPAND DATA WITH ALL POSSIBLE COMBINATIONS:
	.data.long <- .data.long[.line.seq,]
	.navn <- dimnames(.data.long)[[2]] # SAVE NAMES BEFORE RESHAPING
	# REFORMAT DATA TO SET SIDE BY SIDE:
	.data.long <- t(matrix(as.numeric(t(.data.long)), nrow = dim(.data.long)[2]*.nmarkers))
####
####
	if(T){# SOME DATA CHECKING (NOT *REALLY* NECESSARY):
# CHECK LINE NUMBERS:
	.colno.check <- seq(from = 3, by = 5, length = .nmarkers)
	.tmp <- .data.long[,.colno.check, drop = F]-.data.long[,.colno.check[1]]
##	print(.tmp[1:10,])
	if(any(.tmp != 0)) stop("Problem in data preparation!")
#
# CHECK MARKER NUMBERS:
	for (i in 1:.nmarkers){
		.coltmp <- 4 + 5*(i - 1)
		if(any(.data.long[,.coltmp] != i)) stop("Problem in data preparation!")
	}
#
#
	}
	
	# SET CORRECT NAMES AFTER RESHAPING:
	.navn <- as.character(t(outer(paste("l", 1:.nmarkers, sep = ""), .navn, function(x,y) paste(x,y,sep="."))))
	dimnames(.data.long) <- list(NULL, .navn)
#
##
## REORDER COLUMNS AND FIX COLNAMES:
#
###	.indtmp <- as.numeric(outer(7*(0:(.nmarkers - 1)), 1:4,  "+"))
	.indtmp <- as.numeric(outer(5*(0:(.nmarkers - 1)), 1:2,  "+"))
###	.data.long <- .data.long[,c(.indtmp, 5)]
	.data.long <- .data.long[,c(.indtmp, 3)]
	dimnames(.data.long)[[2]][dimnames(.data.long)[[2]] == "l1.ind.unique.line"] <- "ind.unique.line"
#
##
##	COMPUTE HAPLOTYPE NUMBER FOR MARKER COMBINATIONS:
#
	.pos.haplo <- f.pos.to.haplocomb(A = .nalleles, comb = .data.long[,-(dim(.data.long)[2])], fam = "c")
	.haplo.comb <- f.pos.to.haplocomb(pos = .pos.haplo, A = prod(.nalleles), fam = "c")

## return(cbind(.data.long, .pos.haplo, .haplo.comb))

##	f.vis(dimnames(.haplo.comb)[[2]], vis = T)
	dimnames(.haplo.comb)[[2]] <- c("c1", "c2")

#
##
## PREPARE OUTPUT WITH HAPLOTYPE NUMBERS, INDEX OF UNIQUE LINES ()
## AND THE CORRESPONDING FREQUENCIES:
## NOTE: ind.unique.line REFERS TO THE SEQUENCE 1:.nlines.unique () AND NOT THE ORIGINAL LINE NUMBERS
#
##	.ut <- cbind(comb = .haplo.comb, ind.unique.line = .data.long[,"ind.unique.line"])
##	.ut <- cbind(.ut, freq = .data.agg$freq[.ut[,"ind.unique.line"]])


### BOER SJEKKE DENNE?:
##
##
	.ut <- cbind(comb = .haplo.comb, ind.unique.line = .data.long[,"ind.unique.line"], freq = .data.agg$freq[.data.long[,"ind.unique.line"]])
	.ut <- as.data.frame(.ut)
#
	attr(.ut, "alleles") <- .alleles
###	attr(.ut, "rows.with.Mendelian.inconsistency") <- .rows.with.Mendelian.inconsistency
	attr(.ut, "orig.lines") <- .lines # NOTE: .lines CAN BE INDEXED BY .tag AND .tag.unique (AS CHARACTER) OR BY ind.unique.line 
#	
	return(.ut)

}
