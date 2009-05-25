"f.thin"<-
function(data, join, selection)
{
## REMOVES FROM DATA (OR JOINS) SELECTED HAPLOTYPES, OR FIND ZERO (LOW) FREQUENCY 
## ALLELES/HAPLOTYPES AND REMOVE/JOINS THOSE, USING 
## join: SHOULD RARE HAPLOTYPES BE JOINED OR DELETED?
## selection: LOGICAL VECTOR WITH NAMES EQUAL TO THE HAPLOTYPES, TELLING WHAT HAPLOTYPES TO SELECT
##
## IMPORTANT: ASSUMES m1, m2, f1, f2 TO HAVE NUMERIC CODING!
# 
#
	.haplotypes <- attr(data, "haplotypes") #
#
#
	.ind <- which(selection)
	.ind.full <- is.element(data$m1, .ind) & is.element(data$m2, .ind) & is.element(data$f1, .ind) & is.element(data$f2, .ind) #
#
# WARNING: EXCLUDING HAPLOTYPES LIKE THIS SHOULD STILL RESULT IN A DATA SET CORRECTLY SORTED. #
# HOWEVER, IT IS SORTED BY THE ALLELE COMPONENTS OF THE HAPLOTYPES, NOT BY THE ALPHABETIC #
# ORDER OF THE HAPLOTYPES (WELL, THIS HAS BEEN FIXED NOW, BY USING NUMERIC CODES...)#
#
	.ut <- data
	.ut$tmp <- NULL # WARNING: AFTER THINNING THE tmp VARIABLE IS NO LONGER CORRECT, THUS REMOVED
#
#
	if(!join | all(selection)){# NO JOINING OR NOTHING TO JOIN
	.ut <- .ut[.ind.full,] #
#
	attr(.ut, "haplotypes") <- .haplotypes
##	names(.ind0) <- .haplotypes
	attr(.ut, "selected.haplotypes") <- selection
#
	} # END IF NO JOINING
#
#
	else{# JOIN!



	.joincode <- min(which(!selection))
	f.vis("Den nye ambiguity-var her vil ikke stemme for samlet haplotype!\n")
	cat("Tror du virkelig dette vil funke for single locus?? De to amb-variablene?\n")
	for(i in 1:4){
		.ut[!is.element(.ut[,i], .ind),i] <- .joincode	
		}
	
	
	.ut$tmp <- NULL

#
#
	.paste <- paste("#", paste(.ut$m1, .ut$m2, .ut$f1, .ut$f2, sep = "#"), "#", sep="")

	cat("loesriver de nye!\n")


	.ut$freq <- f.groupsum(.ut$freq, .paste) # SUM OVER EACH NEW HAPLOTYPE

	.ut$amb[!.ind.full] <- .paste[!.ind.full]

	.ut <- .ut[!duplicated(.paste),]

	attr(.ut, "haplotypes") <- .haplotypes
##	names(.ind0) <- .haplotypes
	attr(.ut, "selected.haplotypes") <- selection
	cat("Haplotypes", .haplotypes[!selection], "all joined in", .haplotypes[.joincode], "\n")
	
	attr(.ut, "selected.haplotypes")[.joincode] <- T

	if(dim(.ut)[1] != (length(.ind) + 1)^4) stop("Incorrect grid dimension in f.thin\n")
	.tmp <- sort(c(.joincode, .ind))
	if(any(.ut[, 1:4] != expand.grid(.tmp, .tmp, .tmp, .tmp))) stop("Problem in f.thin: Incorrect grid!\n")
#
	}# END if join
#
	return(.ut) 
}
