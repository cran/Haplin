f.sep.data <- function(data, info){
##
## SPLIT GENETIC DATA FROM VARIABLES.
## RESHAPE GENETIC DATA INTO COMPLETE HAPLOTYPE DATA
##
design <- info$model$design
xchrom <- info$model$xchrom
n.vars <- info$filespecs$n.vars
.alleles <- info$haplos$alleles
#
if(n.vars > 0){
	.data.gen <- data[,-(1:n.vars), drop = F] # SELECTS GENETIC DATA
	attr(.data.gen, "alleles") <- .alleles # MAKE SURE TO KEEP ATTRIBUTES
	#
	.data.vars <- data[,1:n.vars, drop = F] # SELECT VARIABLES
	attr(.data.vars, "variables") <- attr(data, "variables")
}else{
	.data.gen <- data # NO CHANGE
	.data.vars <- NULL
}
#
## TEST HWE
.HWE.res <- f.HWE.design(.data.gen, design = design, xchrom = xchrom)
#
## RESHAPE GENETIC DATA
#
## TRIAD + CC.TRIAD DESIGN (NO X-CHROM):
if((design %in% c("triad", "cc.triad")) & !xchrom){
	.data.gen <- f.sort.alleles.new(.data.gen)
}
#
## TRIAD + CC.TRIAD DESIGN (ON THE X-CHROM):
if((design %in% c("triad", "cc.triad")) & xchrom){
	.data.gen <- f.sort.alleles.new(.data.gen, xchrom = T, sex = data[, info$variables$sex])
}
#
## CASE-CONTROL DESIGN:
if(design == "cc"){
	.data.gen <- f.sort.alleles.cc(.data.gen)
}
#
#


return(list(data.gen = .data.gen, data.vars = .data.vars, HWE.res = .HWE.res))





}
