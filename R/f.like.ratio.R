f.like.ratio <- function(res.0, res, data, info){

#
#
## OBTAIN LOG-LIKELIHOOD, VIA DIFFERENT METHODS (SHOULD BE SORTED OUT...):
.anova <- anova(res.0$result, res$result, test = "Chisq")
.df <- .anova$Df[2]



if(F){
	.lratio.test <- anova(res.0$result, res$result, test = "Chisq")[2,"P(>|Chi|)"]
}	
f.vis(.loglike.ratio.full <- .anova$Deviance[2], vis = F)

.loglike.0 <- f.final.loglike(data = data, pred = res.0$pred, info = info)
.loglike <- f.final.loglike(data = data, pred = res$pred, info = info)

f.vis(.loglike.ratio <- 2*(.loglike - .loglike.0), vis = F)

.loglike.ratio <- c(.loglike.ratio, full.loglike = .loglike.ratio.full)
names(.loglike.ratio) <- paste(names(.loglike.ratio), ".ratio", sep = "")

.lratio.test <- 1 - pchisq(.loglike.ratio, df = .df)

.lratio.test <- .lratio.test[1] ## BRUKER BARE DEN FORSTE, DE ANDRE VAR EKSPERIMENTELLE

.lratio.ut <- c(loglike.0 = unname(.loglike.0[1]), loglike = unname(.loglike[1]), df = .df, p.value.overall = unname(.lratio.test))
	
## return(.lratio.test)	

if(F){	
	f.vis(.tmpan <- anova(res.0$result, res$result, test = "Chisq"), vis = T)
	print(anova(res.0$result, res$result, test = "Chisq"), digits = 12)


	f.vis(res.0$result$dev - res$result$dev, vis = T)
}




return(.lratio.ut)


}
