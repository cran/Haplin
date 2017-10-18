f.debug.pvalues <- function(.effects.CI, .pvalues.new, .pvalues.old, .arg.effects, .names.haplo){
##
##
## Debug for calculation of p-values when using Johnson's distribution
##
#
## Create directory for output
dir.create("debug.haplin.pvalues", showWarnings=F)
#
pdf(paste("debug.haplin.pvalues/",paste(.names.haplo,sep="",collapse="-"),".pdf",sep=""))
on.exit(dev.off(),add=T)
par(mfrow=c(3,3))
#
if(length(.arg.effects!=0)){
	.temp.pvalues <- sapply(colnames(.arg.effects),function(x){
		## Johnson's parameters
		.param <- JohnsonFit(log(.arg.effects[,x]), moment="quant")
		#
		## Plots 
		.emp.cum <- ecdf(log(.arg.effects[,x]))
		plot(.emp.cum, main = x)
		.arg.Johnson <- seq(max(min(min(log(.arg.effects[,x]))*1.1,min(log(.arg.effects[,x])*0.9)),-.Machine$double.xmax),min(max(max(log(.arg.effects[,x]))*1.1,max(log(.arg.effects[,x]))*0.9),.Machine$double.xmax),by=0.01)
		.pJohnson <- pJohnson(.arg.Johnson,.param)
		lines(.arg.Johnson,.pJohnson, col="red")
	})
}
#
## Writes estimates, confidence intervals, new and old p-values and their differences to file
write.table(cbind(.effects.CI, p.value.new = as.numeric(.pvalues.new), p.value.old = as.numeric(.pvalues.old), diff.pvalues = (as.numeric(.pvalues.new)-as.numeric(.pvalues.old))), paste("debug.haplin.pvalues/",paste(.names.haplo,sep="",collapse="-"),".txt",sep=""))
#
## Return empty
return(invisible())
#
}