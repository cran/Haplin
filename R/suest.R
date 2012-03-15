suest <- function(reslist){

#
## CHECK FOR ERRORS IN INDIVIDUAL HAPLIN RUNS. REMOVE AND REPORT.
.fjern <- sapply(reslist, function(x) !is.null(attr(x, "error.message")))
.reslist <- reslist[!.fjern]
if(sum(.fjern) > 0) cat('The following estimation results were removed from the list before applying suest \n(due to incomplete estimation results):\n', '"', paste(names(reslist)[.fjern], collapse = '", "'), '"', '\n', sep = "")


.suest <- f.suest(.reslist, debug = F, diag.plot = F)




class(.suest) <- "suest"

return(.suest)






}
