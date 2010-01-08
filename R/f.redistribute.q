f.redistribute <- function(pred, data, design, xchrom){
##
## REDISTRIBUTE OBSERVED FREQUENCIES ACCORDING TO PREDICTED
##
#
##
## MATCH PREDICTED FREQUENCIES TO data:
.n.haplo <- sum(attr(data, "selected.haplotypes"))

if(F){# ERSTATTET AV f.pos.match
	if((design == "triad") & !xchrom){
		.pos <- f.pos.in.grid(A = rep(.n.haplo, 4), comb = as.matrix(data[,c("m1", "m2", "f1", "f2")]))
	}
	if((design == "triad") & xchrom){
		###stop("vet ikke helt her....")
		.pos <- f.pos.in.grid(A = c(rep(.n.haplo, 3), 2), comb = as.matrix(data[,c("m1", "m2", "f2", "sex")]))
		#.pos <- f.pos.in.grid(A = rep(.n.haplo, 4), comb = as.matrix(data[,c("m1", "m2", "f1", "f2")]))
	}
	if(design == "cc"){
		if(xchrom)stop("Not implemented!")
		.pos <- f.pos.in.grid(A = c(rep(.n.haplo, 2), 2), comb = as.matrix(data[,c("c1", "c2", "cc")]))
	}
	if(design == "cc.triad"){
		if(xchrom)stop("Not implemented!")
		.pos <- f.pos.in.grid(A = c(rep(.n.haplo, 4), 2), comb = as.matrix(data[,c("m1", "m2", "f1", "f2", "cc")]))
	}

	.pos.test <- f.pos.match(data = data, design = design, xchrom = xchrom, n.sel.haplos = .n.haplo)
	if(!all.equal(.pos, .pos.test)) stop()
}
.pos <- f.pos.match(data = data, design = design, xchrom = xchrom, n.sel.haplos = .n.haplo)
.pred <- pred[.pos]
	
	
	
###
####
###	MERK: cat("dette kunne vaert gjort en gang for alle...\n")
###



#
##
## RESCALE PREDICTED FREQUENCIES WITHIN EACH TRIAD:
	.predsum <- f.groupsum(X = .pred, INDICES = data$ind)
	.freqsum <- f.groupsum(X = data$freq, INDICES = data$ind) # MERK: .freqsum ER 1 FOR DENNE VARIANTEN
	.pred.redist <- .pred/.predsum * .freqsum
##	.pred.redist <- ifelse(.predsum > 0.00001, .pred/.predsum*f.groupsum(X = data$freq, INDICES = data$ind), 0)
#	
##
## AGGREGATE TRIAD CONTRIBUTIONS OVER HAPLOTYPE COMBINATIONS:
	.pred.redist <- tapply(.pred.redist, .pos, sum)
#	
##
## PREPARE OUTPUT:
	.utfreq <- pred
	.utfreq[] <- 0	
	.utfreq[as.numeric(names(.pred.redist))] <- .pred.redist
#
	.utfreq
}
