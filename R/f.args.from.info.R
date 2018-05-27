f.args.from.info <- function(info){
##
## EXTRACTS VALID haplin ARGUMENTS FROM AN info OBJECT
##
## NOTE: ASSUMES ALL PARAMETERS ARE ON LEVEL 2 IN info, EXCEPT filename
##
#
## EXTRACT STANDARD HAPLIN ARGUMENTS WITH DEFAULTS
.formals <- formals(haplin)
#
## FLATTEN info OBJECT, ONE LEVEL
which.data.par <- which( names( info ) == "data" )
.tmp <- unlist(info[ -which.data.par ], recursive = F)
.tmp <- c( .tmp, data = info[ which.data.par ] )
#
## EXTRACT NAMES FROM ONE LEVEL DOWN
.names <- sapply(info[ -which.data.par ], names)
.names <- c( .names, data = "data" )
.names <- unlist( .names[-1] )
# .names <- c("filename", unlist(.names[-1]))
#
## CHECK NAME LENGTH
if(length(.names) != length(.tmp)) stop( "Problem with args...", call. = FALSE )
#
## USE NAMES FROM ONE LEVEL DOWN
names(.tmp) <- .names
#
## PICK ONLY THOSE THAT ARE ARGUMENTS TO HAPLIN
.tmp <- .tmp[names(.formals)]
#
## FOR SAFETY'S SAKE
.test <- f.check.pars( .tmp, .formals )
#
##
return(.tmp)
}
