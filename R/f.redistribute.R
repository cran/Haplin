f.redistribute <- function(pred, data, pos, freqsum, expand = T){
##
## REDISTRIBUTE OBSERVED FREQUENCIES ACCORDING TO PREDICTED
## NOTE: pred IS OF LENGTH MATCHING THE DESIGN MATRIX, NOT data.
## IT IS TYPICALLY THE RESULT OF THE PREDICTION FROM f.tri.glm,
## I.E. .res$pred
## IF expand = T THE OUTPUT IS OF LENGTH pred. IF NOT, IT 
## FITS THE SIZE OF data
##
#
## MATCH PREDICTED FREQUENCIES TO data:
.pred <- pred[pos]
#
## RESCALE PREDICTED FREQUENCIES WITHIN EACH TRIAD:
.predsum <- f.groupsum(X = .pred, INDICES = data$orig.lines)
.pred.redist <- .pred/.predsum * freqsum
#
if(!expand) return(.pred.redist)
#	
##
## AGGREGATE TRIAD CONTRIBUTIONS OVER HAPLOTYPE COMBINATIONS:
.pred.redist <- tapply(.pred.redist, pos, sum)
#	
##
## PREPARE OUTPUT:
.utfreq <- pred
.utfreq[] <- 0	
.utfreq[as.numeric(names(.pred.redist))] <- .pred.redist
#
return(.utfreq)
}
