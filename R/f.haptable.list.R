f.haptable.list <- function(object){
##
## APPLIES haptable TO EACH ELEMENT OF object
## object IS A list OF haplin OBJECTS
## 
#
#
## CREATE STANDARD haptable FOR EACH ELEMENT
.tab <- lapply(object, haptable)

# this info contains only the information common to all elements in
#   'object' (i.e., either all windows if 'object' is haplinSlide,
#   or all strata if 'object' is haplinStrat)
.info.all <- attr( .tab[[ 1 ]], "info" )
.info.all <- .info.all[ c( "filespecs", "model", "variables", "control" ) ]
class( .info.all ) <- "info"
#
## STACK TABLES (EFFICIENTLY)
.tab <- toDataFrame(.tab)
attr( .tab, "info" ) <- .info.all
#
class(.tab) <- c("haptable", "data.frame")
return(.tab)
}
