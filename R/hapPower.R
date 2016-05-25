hapPower <- function(hapRun.result, alpha = 0.05){
##
## Power simulations
##
## 
#
## Extract information from hapRun object
.hap.run <- hapRun.result
.hapfunc <- attributes(hapRun.result)$hapfunc
#
## Extract p-values from .hap.run
if(.hapfunc == "haplin"){
	.n.sim <- nrow(unique(.hap.run["sim.no"]))
	if(any(is.na(.hap.run["haplos"]))) .hap.run <- .hap.run[-which(is.na(.hap.run["haplos"])),]
	.no.files <- nrow(unique(.hap.run["sim.no"]))
	.hap.run <- cbind(.hap.run["haplos"],.hap.run["reference"],.hap.run["pv.overall"],.hap.run[grep("p.value",names(.hap.run))])
	.reference <- unique(.hap.run[which(.hap.run["reference"]=="ref"),which(colnames(.hap.run)=="haplos")])
	if(length(.reference)> 1) stop("More than one haplotype have been chosen as reference", call. = F)
	.hap.run <- .hap.run[,-which(colnames(.hap.run)=="reference")]
	.pvalues <- lapply(unique(.hap.run$haplos), function(x){ 
		.pvalues <- .hap.run[.hap.run$haplos == x,]
		.pvalues$haplos <- NULL
		.pvalues
	})
	.power <- cbind(haplos=unique(.hap.run$haplos),as.data.frame(t(sapply(.pvalues, function(x) colMeans(x <= alpha, na.rm = TRUE)))))
	.power[.power$haplos==.reference,which(is.na(.power[.power$haplos==.reference,]))] <- "ref"
	if(any(is.na(.hap.run[-which(.hap.run$haplos==.reference),]))) warning("Some p-values were estimated to NA. These were removed from the power calculations", call. = F)
}else if(.hapfunc == "haplinSlide"){
	.n.sim <- length(.hap.run)
	.pvalues <- sapply(.hap.run, function(x){
		if(is.null(attr(x, "Error"))) x$pval.alt.corr
		else NA
	})
	.no.files <- sum(sapply(.hap.run, function(x){
		if(is.null(attr(x, "Error"))) TRUE
		else FALSE
	}))
	if(length(which(!is.na(.pvalues)))!= .no.files)	warning("Some p-values were estimated to NA. These were removed from the power calculations", call. = F)
	.power <- mean(.pvalues <= alpha, na.rm = TRUE)
}else if(.hapfunc == "haplinStrat"){
	.n.sim <- length(.hap.run)
	.error <- sapply(.hap.run, function(x){
		if(is.null(attr(x, "Error"))) TRUE
		else FALSE
	})
	.hap.run[which(!.error)] <- NULL
	.pvalues <- lapply(.hap.run, function(x){
		.pvalues <- x$pval
		names(.pvalues) <- x$test
		.pvalues <- .pvalues[-which(names(.pvalues)=="haplo.freq")]
	}) 
	.pvalues <- do.call("cbind",.pvalues)
	.no.files <- length(.hap.run)
	if(any(is.na(.pvalues))) warning("Some p-values were estimated to NA. These were removed from the power calculations", call. = F)
	.power <- rowMeans(.pvalues <= alpha, na.rm = TRUE)
}	 
#
cat(paste("\nThe power was calculated using ", .no.files, " of ", .n.sim, " files\n", sep = ""))
#
return(.power)
}