print.suest <- function(x, verbose = F, ...){

.suest <- x


cat("\nIndividual score p-values:\n")
print(.suest$pval.obs)
cat("\nMultiple testing corrected p.value, based on minimum of individual values:\n")
print(.suest$pval.obs.corr)
cat("\nNumber of data lines used in each estimation:\n")
print(.suest$lines.account)

if(verbose){

	cat("\n###################\n")
	cat("Stuff used for testing, ignore:\n")
	print(str(.suest[c("bonferroni", "kill", "pval.alt", "pval.alt.corr")]))
}

.test <- cbind("Individual score p-values:", .suest$pval.obs)
#.test <- rbind(.test, .test)



#return(.test)
dimnames(.test) <- list(rep("", dim(.test)[1]), rep("", dim(.test)[2]))


#print(.test, quote = F)

return(invisible(.test))

}
