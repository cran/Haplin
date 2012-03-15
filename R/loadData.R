loadData <- function(database, markers){
##
## LOAD DATA FROM A "HAPLIN" DATABASE
##
#
## GET RELEVANT DATABASE INFO
.file.varnames <- paste(database, "/000varnames.RData", sep = "")
.file.info <- paste(database, "/000info.RData", sep = "")
.file.dumpinfo <- paste(database, "/000dumpinfo.RData", sep = "")
#
## LOAD info
info <- within(list(), load(.file.info))$info
.markers.old <- info$filespecs$markers
if(missing(markers) || identical(markers, "ALL")) markers <- .markers.old
if(any(!is.element(markers, .markers.old))) stop("Invalid markers selected", call. = F)
#
.n.vars <- info$filespecs$n.vars
.fam <- info$model$fam
#
## LOAD VARIABLE NAMES varnames
varnames <- within(list(), load(.file.varnames))$varnames
#
## SELECT COLUMNS TO USE
.cols <- f.sel.markers(n.vars = .n.vars, markers = markers, family = .fam, split = T, ncols = length(varnames))
#
## LOAD DATA FROM DATABASE
.data <- f.load(database = database, n = length(varnames), cols = .cols, varnames = varnames)
#
##
return(.data)
}

