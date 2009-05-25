f.redistribute <- function(pred, data, design){
##
## REDISTRIBUTE OBSERVED FREQUENCIES ACCORDING TO PREDICTED
##
#
##
## MATCH PREDICTED FREQUENCIES TO data:
	.n.haplo <- sum(attr(data, "selected.haplotypes"))

if(design == "triad"){
	.pos <- f.pos.in.grid(A = rep(.n.haplo, 4), comb = as.matrix(data[,c("m1", "m2", "f1", "f2")]))
}
if(design == "cc"){
	.pos <- f.pos.in.grid(A = c(rep(.n.haplo, 2), 2), comb = as.matrix(data[,c("c1", "c2", "cc")]))
}
if(design == "cc.triad"){
	.pos <- f.pos.in.grid(A = c(rep(.n.haplo, 4), 2), comb = as.matrix(data[,c("m1", "m2", "f1", "f2", "cc")]))
}
	.pred <- pred[.pos]
	
	
	
###
####
###	MERK: cat("dette kunne vaert gjort en gang for alle...\n")
###




if(F){# BARE EN TEST, IKKE NODVENDIG:
	.dim1 <- round(length(pred)^(1/4))
	f.vis(.vekt <- .dim1^(0:3))	
	f.vis(.pos.old <- (data$m1 - 1) + (data$m2 - 1)*.vekt[2] + (data$f1 - 1)*.vekt[3] + (data$f2 - 1)*.vekt[4] + 1)
	if(any(.pos != .pos.old)) stop()
##	cat("bare til test!\n")
}



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
