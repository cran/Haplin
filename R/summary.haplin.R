summary.haplin <- function(object, reference, ...){ #


#
## info OBJECT
.info <- object$info
.poo <- .info$model$poo
#
## Previously used reference method
.reference.method <- .info$haplos$reference.method
.selected.haplotypes <- .info$haplos$selected.haplotypes
#
## USE APPROPRIATE REFERENCE
if(!missing(reference)){
	if(reference == "reciprocal" & .poo){
		warning(paste('Can only (for the time being) use reference = "ref.cat" or "population" when poo == TRUE. Has been changed to ', .reference.method, sep = ""), call. = F)
	} else
	if(reference %in% c("reciprocal", "population", "ref.cat")){	
		.reference.method <- reference
	} else if (is.numeric(reference)){
		cat("\nWARNING: REFERENCE CATEGORY CAN ONLY BE SET IN FIRST RUN OF HAPLIN!\nFOR summary AND plot METHODS ONLY REFERENCE METHOD CAN BE CHOSEN, NOT CATEGORY\n\n")
	} else stop("Invalid reference choice!", call. = F)
}
#
## CHECK THAT ONLY REFCAT IS USED WHEN ONLY TWO HAPLOTYPES/ALLELES	
if(sum(.selected.haplotypes) == 2 & .reference.method != "ref.cat"){
	cat("\nNOTE: ONLY SINGLE REFERENCE CATEGORY METHOD ALLOWED FOR TWO HAPLOTYPES/ALLELES!\n (reference has been set to", object$result$ref.cat, ")\n")
	.reference.method <- "ref.cat"
	}
#
#
.summ.res <- summary.tri.glm(object$result, reference.method = .reference.method, info = .info, ...) #

.ut <- list(summary.tri.glm = .summ.res, info = .info, alleles = .info$haplos$alleles, selected.haplotypes = .selected.haplotypes, HWE.res = .info$check$HWE.res, ntri.seq = .info$data$ntri.seq, loglike = object$loglike, score = object$score)
class(.ut) <- "summary.haplin"

return(.ut)
}
