lineByLine <- function(infile, outfile, linefunc = identity , choose.lines = NULL, choose.columns = NULL, ask = TRUE, blank.lines.skip = TRUE, verbose = TRUE, ...){
##
## Reads file line by line, modifies each line using the argument linefunc, and then writes to outfile.
##
## By default, linfunk returns its argument.
## Default sets ask = TRUE. No file is overwritten unless specified by user.
## choose.lines and choose.columns enable the selection of certain lines and/or columns of the infile. 
## Both are set to NULL as default, which means that all lines and/or columns are selected. If not set to NULL, must be numeric vectors with values > 0.
## blank.lines.skip ignores blank lines in the input if TRUE.
## If TRUE, verbose displays the line number for each iteration, in addition to output from linefunc.
##
#
## Error in lineByLine if outfile == infile
if(file.exists(outfile) & outfile == infile) stop("\"outfile\" is equal to \"infile\"", call. = F)
#
## Open infile for reading and writing
.infile <- file(description = infile, open = "r+")
#
## Make sure the connection is closed when function exits
on.exit(close(.infile))
#
## If outfile exists & ask == TRUE, query user
if(file.exists(outfile) & ask){
	.answer <- readline(paste('Overwrite ', outfile, '? (y/n)', sep = ""))
	if(.answer != "y"){
		cat("Stopped without overwriting file\n")
		return(invisible())
	}
}
unlink(outfile)
#
## (Re)open outfile for writing
.outfile <- file(description = outfile, open = "w")
#
## Make sure the connections are closed when function exits
on.exit(close(.outfile), add = TRUE)
#
## Error if choose.lines has duplicated values
if(any(duplicated(choose.lines))) stop("\"choose.lines\" contains duplicated values", call. = F)
#
## Loop over lines
.i <- 0
#
.k <- 0 # Equals line number if line does not have enough columns 
#
repeat{
	.i <- .i + 1
	#
	## Break off at a given number of lines
	if(!is.null(choose.lines)){
		if(.i > max(choose.lines) + 0.1) break
	}
	#
	## Read a single line. NOTE: Reads numeric as character
	.line <- readLines(.infile, n = 1) 
	#
	## Break off at end of file
	if(length(.line) == 0) break 
	#
	## Skip blank lines if blank.lines.skip == TRUE. Otherwise a warning is given
	if(!nzchar(.line) & blank.lines.skip){
		.i <- .i-1
		next
	}
	if(!nzchar(.line) & !blank.lines.skip){
		warning(paste("Line ", .i, " is empty", sep = ""), call. = F)
	}
	#
	## Choose lines 
	if(!is.null(choose.lines)){
		if(!is.numeric(choose.lines)) stop("\"choose.lines\" must be numeric", call. = F)
		if(sum(choose.lines < 0) != 0) stop("Invalid line number(s)", call. = F)
		if(!is.element(.i, choose.lines)) next
	}	
	#
	## Split line into elements	
	.line <- strsplit(.line, split = " ", fixed = T)[[1]]
	#
	## Choose columns, allowing reordering
	if(!is.null(choose.columns)){
		if(!is.numeric(choose.columns)) stop("\"choose.columns\" must be numeric", call. = F)
		if(sum(choose.columns > 0) != length(choose.columns)) stop("Invalid column number(s)", call. = F)
		if(sum(choose.columns <= length(.line)) != length(choose.columns))	.k = .i
		.line <- .line[choose.columns[which(choose.columns <= length(.line) & choose.columns > 0)]]
	}	
	#
	## Display line number (and output from linefunc) 
 	if(verbose) cat(.i, " --- ", sep = "")
	#
	## Convert line
	if("verbose" %in% names(formals(linefunc))){
		.line <- linefunc(x = .line, verbose = verbose, ...)
	}else{
		.line <- linefunc(x = .line, ...)
	}
	#
	## Newline
	if(verbose) cat("\n")
	#
	## Display invalid column numbers
	if(.k == .i){
		if(verbose) cat("Invalid column number(s).\n")
	}
	#
	## Paste elements, separated by space
	.line <- paste(.line, collapse = " ")
	#
	## Skip blank lines in order to delete the requested rows/lines
	if(nzchar(.line) == FALSE) next
	#
	## Write to new file
	writeLines(.line, .outfile) 
}# end repeat
#
## Display number of lines read and converted
if(is.null(choose.lines)){ 
	cat("Read and converted information for", .i-1, "line(s).\n", sep = " ")
}else if(!is.null(choose.lines) & (sum(choose.lines <= .i-1) == length(choose.lines))){
	cat("Read and converted information for", length(choose.lines), "line(s).\n", sep = " ")
}else{
	warning(paste("Invalid line number(s). Read and converted information for", sum(choose.lines <= .i-1), "line(s).\n", sep = " "), call. = F)
}	
#
## Warning if choose.column contains invalid column number(s)
if(.k != 0) warning("\"choose.columns\" contains invalid column number(s).", call. = F)
#
## Return number of lines read
return(invisible(.i - 1))
}

