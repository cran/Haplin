f.pos.match <- function(data, design, xchrom, n.sel.haplos){
##
## FOR data, COMPUTES THE POSITION OF EACH LINE IN A COMPLETE
## GRID AS USED IN THE GLM
#
##
## PREPARE:
###.n.sel.haplos <- sum(info$haplos$selected.haplotypes)
###.design <- info$model$design
###.xchrom <- info$model$xchrom
.n.sel.haplos <- n.sel.haplos
.design <- design
.xchrom <- xchrom
#
##
## MATCH PREDICTED PROBABILITIES TO ORIGINAL DATA:
if((.design == "triad") & !.xchrom){
	.pos <- f.pos.in.grid(A = rep(.n.sel.haplos, 4), comb = as.matrix(data[,c("m1", "m2", "f1", "f2")]))
}
if((.design == "triad") & .xchrom){
	.pos <- f.pos.in.grid(A = c(rep(.n.sel.haplos, 3), 2), comb = as.matrix(data[,c("m1", "m2", "f2", "sex")]))
}
if(.design == "cc"){
	if(.xchrom)stop("Not implemented")
	.pos <- f.pos.in.grid(A = c(rep(.n.sel.haplos, 2), 2), comb = as.matrix(data[,c("c1", "c2", "cc")]))
}
if(.design == "cc.triad"){
	if(.xchrom)stop("Not implemented")
	.pos <- f.pos.in.grid(A = c(rep(.n.sel.haplos, 4), 2), comb = as.matrix(data[,c("m1", "m2", "f1", "f2", "cc")]))
}
return(.pos)
}
