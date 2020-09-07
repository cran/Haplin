haptable.gxe <- function(object){
##
## Extract part of the gxe test result 
## and put in a format suitable to be 
## joined with a haptable
##
.tab <- object$gxe.test
.gxe.test <- .tab$gxe.test
.tab <- .tab[, "pval", drop = F] # for the time being, pick only pval. But others can be added directly here
#
## Colnames
.colnames <- outer(.gxe.test, colnames(.tab), paste, sep = "_")
.colnames <- as.character(t(.colnames))
.colnames <- paste0("gxe_", .colnames)
#
## Contents
.values <- matrix(t(as.matrix(.tab)), nrow = 1)
#
## Return as data frame
.values <- as.dframe(.values)
colnames(.values) <- .colnames
#
##
return(.values)
}
