f.sep.data <- function(data, info){

design <- info$model$design
xchrom <- info$model$xchrom
n.vars <- info$filespecs$n.vars
.alleles <- info$haplos$alleles

.data <- data



# TRIAD DESIGN (NO X-CHROM):
	if(design == "triad" & !xchrom){
		if(n.vars == 0){
			.data.gen <- .data # NO CHANGE
		}
		if(n.vars > 0){
			.data.gen <- .data[,-(1:n.vars), drop = F] # SELECTS GENETIC DATA
			attr(.data.gen, "alleles") <- .alleles # MAKE SURE TO KEEP ATTRIBUTES
		}
		## TEST FOR HWE
		.HWE.res <- f.HWE.design(.data.gen, design = design)
		## RESHAPE DATA
		.data.gen <- f.sort.alleles.new(.data.gen)
	}# END TRIADS
#
## TRIAD DESIGN (ON THE X-CHROM):
	if(design == "triad" & xchrom){
		if(n.vars == 0){
			stop("Need a variable to specify sex!\n")
			#.data.gen <- .data # NO CHANGE
		}
		if(n.vars > 0){
			.data.gen <- .data[,-(1:n.vars), drop = F] # SELECTS GENETIC DATA
			attr(.data.gen, "alleles") <- .alleles # MAKE SURE TO KEEP ATTRIBUTES
		}
		## TEST FOR HWE
		.HWE.res <- f.HWE.design(.data.gen, design = design)
		## RESHAPE DATA
		.data.gen <- f.sort.alleles.new(.data.gen, xchrom = T, sex = .data[,1])
		cat("ADVARSEL: KJOENNSVARIABEL SATT TIL 1!\n")
	}# END TRIADS
#
## CASE-CONTROL DESIGN:
	if(design == "cc"){
		if(xchrom) stop("Not yet implemented")
		if(n.vars == 0){
			stop("Problem: case-control data need at least a variable indicating case-control status! (n.vars must be greater than 0)")
		}
		if(n.vars > 0){
			.data.gen <- .data[,-(1:n.vars), drop = F] # SELECTS GENETIC DATA
			attr(.data.gen, "alleles") <- .alleles # MAKE SURE TO KEEP ATTRIBUTES
			## TEST FOR HWE
			.HWE.res <- f.HWE.design(.data.gen, design = design)
			## RESHAPE DATA
			.data.gen <- f.sort.alleles.cc(.data.gen)
			.data.vars <- .data[,1:n.vars, drop = F] # SELECT VARIABLES
			attr(.data.vars, "variables") <- attr(.data, "variables")
		}
	}# END CASE-CONTROL
#
## COMBINED TRIAD AND CASE-CONTROL:
	if(design == "cc.triad"){
		if(xchrom) stop("Not yet implemented")
		if(n.vars == 0){
			stop("Problem: case-control data need at least a variable indicating case-control status! (n.vars must be greater than 0)")
		}
		if(n.vars > 0){
			.data.gen <- .data[,-(1:n.vars), drop = F] # SELECTS GENETIC DATA
			attr(.data.gen, "alleles") <- .alleles # MAKE SURE TO KEEP ATTRIBUTES
			## TEST FOR HWE
			.HWE.res <- f.HWE.design(.data.gen, design = design)
			## RESHAPE DATA
			.data.gen <- f.sort.alleles.new(.data.gen)
			.data.vars <- .data[,1:n.vars, drop = F] # SELECT VARIABLES
			attr(.data.vars, "variables") <- attr(.data, "variables")
		}
	}# END CASE-CONTROL-TRIADS
#
#

if(!exists(".data.vars")) .data.vars <- NULL

return(list(data.gen = .data.gen, data.vars = .data.vars, HWE.res = .HWE.res))





}
