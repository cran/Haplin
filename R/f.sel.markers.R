f.sel.markers <- function(n.vars, markers, family, split, all.marker.names, ncols){
##
## COMPUTES THE COLUMN NUMBERS OF MARKERS TO BE EXTRACTED FROM A
## FILE IN SEPARATE (OR NOT SEPARATE) COLUMN FORMATS
## NOTE: n.vars MUST BE ADDED, IF NECESSARY
##
#
## SET NUMBER OF COLUMNS PER FAMILY (OR PER CHILD)
if(missing(ncols)){
	ncols <- length(all.marker.names)
}
if(family == "mfc"){
	if(split) .t <- 6
	else .t <- 3
}else if(family == "c"){
	if(split) .t <- 2
	else .t <- 1
} else stop("Problem with family argument!", call. = F)
.test <- (ncols - n.vars) / .t 
.err <- F
#
if(!is.numeric(markers)){
  if(length(markers) == 1 & markers == "ALL"){
	  if(round(.test) != .test) .err <- T
  	.markers <- 1:.test
  	.sel <- 1:ncols
  } else {
    # 'markers' is a character vector with names of markers
    .sel <- grep(
      pattern = paste0(markers, "_", collapse="|"),
      x = all.marker.names,
      value = FALSE
    )
  }
}else{
	## COMPUTE COLUMN NUMBERS OF MARKERS SPECIFIED IN markers ARGUMENT 
	.sel <- c(seq(length.out = n.vars), n.vars + as.numeric(t(outer((markers-1)*.t, 0:(.t - 1), "+")) + 1))
	if(max(.sel) > ncols) .err <- T
	.markers <- markers
}
if(.err){
	stop('It appears that the number of columns in the data file is wrong! \nCheck the datafile, the "n.vars", "sep" and "markers" arguments, etc....', call. = F)
}
#
##
attr(.sel, "markers") <- .markers
attr(.sel, "nloci") <- length(.markers)
return(.sel)
}
