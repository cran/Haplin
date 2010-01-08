f.HWE.design <- function(data, design, xchrom){
##
## TESTS EACH LOCUS FOR HWE USING A STANDARD CHI-SQUARED.
## INTENDED ONLY FOR THE GENETIC PART OF THE DATA, WITH STANDARD TRIAD
## STRUCTURE
##

.alleles <- attr(data, "alleles")
.nloci <- length(.alleles)

if(!is.numeric(data)) stop("Data format wrong when testing HWE!")

.ut <- vector(.nloci, mode = "list")

for(i in 1:.nloci){
	## FOR EACH LOCUS, STACK MOTHER, FATHER, CHILD AND COMPUTE JOINTLY (IF TRIADS)
	if(design %in% c("triad", "cc.triad")){
		.data.tmp <- data[, ((i-1)*6 + 1):((i-1)*6 + 6)]
		if(!xchrom){
			.data.tmp <- rbind(.data.tmp[,1:2], .data.tmp[,3:4], .data.tmp[,5:6])
		}else{
			.data.tmp <- .data.tmp[,1:2] # USES ONLY MOTHERS
		}
	}
	if(design == "cc"){
		.data.tmp <- data[, ((i-1)*2 + 1):((i-1)*2 + 2)]
	}
	.HWE.test <- f.HWE(.data.tmp, quiet.warning = T)
	.ut[[i]] <- .HWE.test
}


return(.ut)

}
