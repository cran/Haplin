"print.summary.tri.glm"<-
function(x, digits = 2, haplos, ...)
{
# PRINTS THE RESULT FROM summary.tri.glm
#
#### PREPARE: #####################################
	.n.all <- x$n.all
	.effects <- x$effects	#
	design <- x$design
	if(missing(haplos)) haplos <- paste("h", 1:.n.all, sep = "")
	.sel.sex <- x$info$variables$sel.sex
	.comb.sex <- x$info$variables$comb.sex
	.response <- x$info$haplos$response
	.xchrom <- x$info$model$xchrom
#
#### PRINT GENERAL INFORMATION: ####################
	cat("\nDate of call:\n")
	cat(x$date, "\n")	
	cat("\nNumber of triads: ", round(x$n.tri), "\n", sep = "")
	cat("\nNumber of haplotypes: ", round(x$n.all), "\n", sep = "")	#
#
# 
#### PRINT ALLELE FREQUENCIES: #####################
if(!x$conf.int){
	cat("\nHaplotype frequencies (%):\n")
	.printout <- 100*round(.effects[1:.n.all, 1], digits)
	names(.printout) <- haplos
	print(.printout, quote = F)	#
}else{
	.printout <- format(100*.effects[1:.n.all, 1:3], digits = digits + 0, width = digits + 2 )
	dimnames(.printout)[[2]][1] <- "Frequency(%)"
	.printout <- cbind(Haplotype = haplos, .printout)
	dimnames(.printout)[[1]][] <- ""		
	cat("\nHaplotype frequencies with ", 100 * x$level, "% confidence intervals:\n", sep = "")
	print(.printout, quote = F)	#
} #
#
#### PRINT RELATIVE RISKS: #######################
.printout <- formatC(.effects[ - c(1:.n.all), , drop = F], flag = "-", digits = digits, width = max(digits + 2, 10) )
#
if(x$reference.method == "ref.cat"){
	## REMOVE PRINTOUT FOR REFERENCE CATEGORY	
	if(x$maternal){
		if(.n.all == 2)
			.ind.strikeout <- x$ref.cat + c(0, 2, 4, 6)
		else
			.ind.strikeout <- x$ref.cat + .n.all * c(0,2)
	}
	else { 
		if(.n.all == 2)
			.ind.strikeout <- x$ref.cat + c(0, 2)
		else
			.ind.strikeout <- x$ref.cat
	}
} else .ind.strikeout <- NULL
#
#
if(!x$conf.int){
	cat("\nSingle- and double dose effects (Relativ Risk):\n")
	.printout[.ind.strikeout,] <- "REF"
}
else {
	cat("\nSingle- and double dose effects (Relativ Risk) with ", 100 * x$level, "% confidence intervals:\n", sep = "")
	if(design == "triad") for (i in seq(along = .ind.strikeout)) .printout[.ind.strikeout[i],] <- c("REF", "", "", "")
}
#
## PRINT REFERENCE METHOD/CATEGORY
cat("Reference method: ", x$reference.method, "\n", sep = "")
if(x$reference.method == "ref.cat") cat("Reference category: ", x$ref.cat, "\n", sep = "")
if(.xchrom){
	cat("Assuming X-chromosome data:\n")
	if(.comb.sex == "males"){
		cat("   |Estimating risk from males only\n")
	}
	if(.comb.sex == "females"){
		cat("   |Estimating risk from females only\n")
	}
	if(.comb.sex == "single"){
		cat(paste('   |With comb.sex = "single", males are assumed to have\n   |a single dose.\n', sep = ""))
	}
	if(.comb.sex == "double"){
		cat('   |With comb.sex = "double", males are assumed to have\n   |a double dose.\n')
	}
}



cat("\n")

if(!x$conf.int){
	dimnames(.printout)[[2]] <- c("Relative Risk")
}
else {
	dimnames(.printout)[[2]] <- c("Relative Risk", "Lower CI", "Upper CI", "P-value")
}
	.printout <- cbind(Haplotype = haplos, .printout)

.flett <- T

if(!.flett){
	## GAMMEL VARIANT, LENGE SIDEN OPPDATERING
	cat("Child single dose:\n")
	print(.printout[1:.n.all,], quote = F)
	cat("\n")
	
	cat("Child double dose:\n")
	print(.printout[(.n.all+1):(2*.n.all),], quote = F)
	cat("\n")
	
	if(x$maternal){
		cat("Maternal single dose:\n")
		print(.printout[(2*.n.all + 1):(3*.n.all),], quote = F)
		cat("\n")
		
		cat("Maternal double dose:\n")
		print(.printout[(3*.n.all + 1):(4*.n.all),], quote = F)
		cat("\n")
	} # END IF MATERNAL
} # END IF NOT FLETT
#
else{
	#
	## RE-ARRANGE SEQUENCE, JUXTAPOSE SINGLE- AND DOUBLE DOSE
	.ind.child <- as.numeric(t(cbind(matrix(1:(2*.n.all), ncol = 2), NA)))
	.printout.child <- cbind(Haplotype = .printout[.ind.child, 1], Dose = c("S  ", "D  ", ""), .printout[.ind.child,-1, drop = F])
	.printout.child[is.na(.printout.child)] <- ""
	dimnames(.printout.child)[[1]][] <- ""
	if(!is.null(.sel.sex) && .sel.sex == 1){
		## REMOVE SINGLE OR DOUBLE DOSE FOR BOYS IN CASE X-CHROM AND SELECTED ONLY BOYS
		.kill <- grep("D", .printout.child[, "Dose"])
		.printout.child <- .printout.child[-.kill, , drop = F]
	}
	#
	## PRINT
	cat("----Child haplotypes----\n")
	print(.printout.child, quote = F)
	#
	##	
	if(x$maternal){
		#
		## RE-ARRANGE SEQUENCE, JUXTAPOSE SINGLE- AND DOUBLE DOSE
		.ind.maternal <- as.numeric(t(cbind(matrix((2*.n.all + 1):(4*.n.all), ncol = 2), NA)))
		.printout.maternal <- cbind(Haplotype = .printout[.ind.maternal, 1], Dose = c("S  ", "D  ", ""), .printout[.ind.maternal,-1, drop = F])
		.printout.maternal[is.na(.printout.maternal)] <- ""
		dimnames(.printout.maternal)[[1]][] <- ""
		#
		## PRINT
		cat("----Maternal haplotypes----\n")
		print(.printout.maternal, quote = F)
	} # END IF MATERNAL
}# END ELSE FLETT
#
#### END: #################################
invisible(x)
}
