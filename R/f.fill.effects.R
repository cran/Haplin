f.fill.effects <- function(resmat, info){
## FILL IN "MISSING" COLUMNS (NON-ESTIMATED EFFECT 
## PARAMETERS) WITH ZEROS
##
#
## PREPARE
.n.sel.haplos <- sum(info$haplos$selected.haplotypes)
.maternal <- info$model$maternal
.resp <- info$haplos$response
#
## BUILDING BLOCKS FOR EFFECT NAMES
.mf <- paste("mf", 1:.n.sel.haplos, sep = "")
.c <- paste("c", 1:.n.sel.haplos, sep = "")
.cdd <- paste("cdd", 1:.n.sel.haplos, sep = "")
.m <- paste("m", 1:.n.sel.haplos, sep = "")
.mdd <- paste("mdd", 1:.n.sel.haplos, sep = "")
#
## SET UP RELEVANT EFFECT NAMES VECTOR:
if(!.maternal){
	if(.resp == "free"){
		.navn <- c(.mf, .c, .cdd)
	}
	if(.resp == "mult"){
#		.navn <- c(.mf, .c)
		.navn <- c(.mf, .c, .cdd)
	}
}
if(.maternal){
	if(.resp == "free"){
		.navn <- c(.mf, .c, .cdd, .m, .mdd)
	}
	if(.resp == "mult"){
#		.navn <- c(.mf, .c, .m)
		.navn <- c(.mf, .c, .cdd, .m, .mdd)
	}
}
#
## CHECK FOR INCORRECT NAMES
.resnavn <- dimnames(resmat)[[2]]
if(any(!is.element(.resnavn, .navn))) stop("Problem with effect matrix")
#
## PAD WITH ZEROS, BY CONVERTING TO LIST AND BACK AGAIN
.ut <- f.matrix.to.list(resmat)
names(.ut) <- .resnavn
.ut <- .ut[.navn]
names(.ut) <- .navn
.ut[!is.element(.navn, .resnavn)] <- 0
.ut <- do.call("cbind", .ut)


f.vis(head(resmat), vis = F)
f.vis(head(.ut), vis = F)

return(.ut)

}
