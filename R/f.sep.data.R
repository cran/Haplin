f.sep.data <- function( data, info ){
##
## RESHAPE GENETIC DATA INTO COMPLETE HAPLOTYPE DATA
##
design <- info$model$design
xchrom <- info$model$xchrom
n.vars <- info$filespecs$n.vars
alleles <- info$haplos$alleles

nloci <- length( alleles ) 
# checking whether any locus contain only one allele
# (this can happen if we omit some data, e.g., due to missingness);
# if so, re-code this data
for(i in 1:nloci){
	alleles.count <- length( alleles[[ i ]] )
	if( alleles.count == 1 ){
		if(design %in% c("triad", "cc.triad")){
			chosen.cols <- ((i-1)*6 + 1):((i-1)*6 + 6)
		}
		if(design == "cc"){
			chosen.cols <- ((i-1)*2 + 1):((i-1)*2 + 2)
		}
		data$gen.data[ ,chosen.cols ] <- 1
	}
}

## TEST FOR HWE
if(!xchrom){
	.HWE.res <- f.HWE.design( data, design = design )
} else {
	.HWE.res <- f.HWE.design( data, design = design, sex = data$cov.data[, info$variables$sex] )
}

## RESHAPE GENETIC DATA
if( xchrom ){
	.data.gen <- f.sort.alleles.new( data, design = design, xchrom = T, sex = data$cov.data[, info$variables$sex] )
} else {
	.data.gen <- f.sort.alleles.new( data, design = design )
}

return( list(data.gen = .data.gen, HWE.res = .HWE.res) )
}
