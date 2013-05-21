convertPed <- function(ped.infile, map.infile, ped.outfile, map.outfile, create.unique.id = FALSE, convert, snp.select = NULL, choose.lines = NULL, ask = TRUE, blank.lines.skip = TRUE, verbose = TRUE){
##
## Recodes a ped file with several options: 
## 	- Possibly creating unique individual IDs
## 	- Possibly converting the SNP alleles between letters (A,C,G,T) and numbers (1,2,3,4) 
## 	- Possibly making a selection of SNPs.
##
## ped.infile is a character string giving the name and path of the standard ped file to be modified.
## map.infile is a character string giving the name and path of the to-be-modified standard map file.
## ped.outfile is a character string of the name and path of the converted ped file.
## map.outfile is a character string giving the name and path of the modified map file.
## snp.select is a numeric vector of SNP numbers or a character vector of the SNP identifiers (RS codes). Default is NULL, which means that all SNPs are selected. 
## ask is a logical variable. If "TRUE", convertPed will ask before overwriting an already existing 'outfile'.
## choose.lines enables the selection of certain lines of the infile. 
## Default is NULL, which means that all lines are selected. If not set to NULL, must be numeric vectors with values > 0.
## blank.lines.skip is a logical varible. If "TRUE", convertPed ignores blank lines in ped.infile and map.infile.
## create.unique.id. If "TRUE", the function creates a unique individual ID.
## convert. The option "ACGT_to_1234" recodes the SNP alleles from A,C,G,T to 1,2,3,4, whereas "1234_to_ACGT" converts from 1,2,3,4 to A,C,G,T. If "no_recode", no conversion occurs.
## verbose. Default is "TRUE", which means that the line number is displayed for each iteration, i.e. each line read and modified, in addition to the first ten columns of the converted line.
##
#
## Error if map.outfile is equal to map.infile
if(file.exists(map.outfile) & map.outfile == map.infile) stop("\"map.outfile\" is equal to \"map.infile\"", call. = F)
# 
## Read map file
.map.infile = read.table(map.infile, header = TRUE, stringsAsFactors = FALSE, blank.lines.skip = blank.lines.skip)
#
## Checking if the map and ped file have the same number of SNPs
.ped.infile <- file(description = ped.infile, open = "r+")
.line <- readLines(.ped.infile, n = 1)
.columns.ped.infile <- length(strsplit(.line, split = " ", fixed = T)[[1]])
close(.ped.infile)
if(.columns.ped.infile != 2*nrow(.map.infile)+6) stop("The map and ped file do not have the same number of SNPs", call. = F)
#
## If snp.select is set to NULL (default), all SNPs are selected
if(is.null(snp.select)) snp.select <- .map.infile[,2]
#
## Find the column positions in the ped file corresponding to the vector snp.select using the function snpPos
.snp.select.num <- snpPos(snp.select = snp.select, map.file = map.infile, blank.lines.skip = blank.lines.skip)
#
## Error if choose.lines has duplicated values
if(any(duplicated(choose.lines))) stop("\"choose.lines\" contains duplicated values", call. = F)
#
## Convert data frame containing map data
if(is.numeric(snp.select)) snp.select <- .map.infile[snp.select,2]
.converted.map <- .map.infile[match(snp.select, .map.infile[,2]),] 
#
## Create new map file if map.outfile does not exist. Overwrite existing file if requested.
if(file.exists(map.outfile) & ask){
	.answer <- readline(paste('Overwrite ', map.outfile, '? (y/n)', sep = ""))
	if(.answer != "y"){
		cat("Stopped without overwriting file\n")
		return(invisible())
	}
}
write.table(.converted.map, map.outfile, row.names = F, col.names = T, quote = F)
#
## Display number of SNPs in map file
cat(paste("Created map file for", length(snp.select), "SNPs \n", sep = " "))
# 
## Line-by-line conversion of ped-file. Message to user
cat("Starting line-by-line conversion of ped file \n")
#
## Convert ped file using the function lineByLine
lineByLine(infile = ped.infile, outfile = ped.outfile, linefunc = lineConvert, choose.lines = choose.lines, ask = ask, blank.lines.skip = blank.lines.skip, verbose = verbose, create.unique.id = create.unique.id, convert = convert, snp.select.num = .snp.select.num)
#
## Return empty
return(invisible())
#
}
