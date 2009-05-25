f.final.loglike <- function(data, pred, info){
##
## COMPUTES THE MAXIMUM LOG-LIKELIHOOD (UP TO A CONSTANT) FOR THE FINAL RESULT
## NOTE: TAKES INTO ACCOUNT MISSING INFORMATION, I.E. IS not THE FULL
## LIKELIHOOD IN EM BUT RATHER THE CORRECT MAXIMUM LIKELIHOOD FROM THE OBSERVED
## DATA
##
## data IS A DATA FRAME WITH A TRIAD INDICATOR SHOWING ALL POSSIBLE HAPLOTYPE
## COMBINATIONS FOR THAT TRIAD, pred ARE THE PREDICTED FREQUENCIES IN THE 
## MAXIMIZED FULL LIKELIHOOD, FOR ALL HAPLOTYPE COMBINATIONS IN STANDARD
## ORDERING
##
##
#
##
## PREPARE:
.n.sel.haplos <- sum(info$haplos$selected.haplotypes)
.design <- info$model$design
#
##
## STANDARDIZE TO PROBABILITIES, AS IN A MULTINOMIAL:
.prob <- pred/sum(pred)
#
##
## MATCH PREDICTED PROBABILITIES TO ORIGINAL DATA:

if(.design == "triad"){
	.pos <- f.pos.in.grid(A = rep(.n.sel.haplos, 4), comb = as.matrix(data[,c("m1", "m2", "f1", "f2")]))
}
if(.design == "cc"){
	.pos <- f.pos.in.grid(A = c(rep(.n.sel.haplos, 2), 2), comb = as.matrix(data[,c("c1", "c2", "cc")]))
}
if(.design == "cc.triad"){
	.pos <- f.pos.in.grid(A = c(rep(.n.sel.haplos, 4), 2), comb = as.matrix(data[,c("m1", "m2", "f1", "f2", "cc")]))
}

.prob <- .prob[.pos]


if(T){

# TESTER LITT HER... DETTE BURDE SVARE TIL DEN FULLE?

	.probsum <- f.groupsum(.prob, data$ind)
	.probnorm <- .prob/.probsum

	.test.loglike <- sum(.probnorm * log(.prob))

}


#
##
## SUM PREDICTED PROBABILITIES OVER AMBIGUITIES FOR EACH TRIAD:
	.prob <- tapply(.prob, data$ind, sum)
#
##
## COMPUTE LOG-LIKELIHOOD:
	.loglike <- sum(log(.prob)) # NOTE: FREQUENCIES ARE 1
#	
##
## FINISH:
return(c(loglike = .loglike, test.loglike = .test.loglike))

}
