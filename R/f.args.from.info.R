f.args.from.info <- function(info){
##
## EXTRACTS VALID haplin ARGUMENTS FROM AN info OBJECT
##

.formals <- formals(haplin)
#
## NOTE: ASSUMES EVERYTHING IS ON LEVEL 2, EXCEPT filename

.tmp <- unlist(info, recursive = F)

.names <- sapply(info, names)
.names <- c("filename", unlist(.names[-1]))

if(length(.names) != length(.tmp)) stop()

names(.tmp) <- .names

.tmp <- .tmp[names(.formals)]

.test <- f.check.pars(.tmp, .formals)

return(.tmp)

}
