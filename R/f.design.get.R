f.design.get <- function(n.all = 2, design = "triad", xchrom = F, maternal = F, poo = F, hwe = T, comb.sex = NULL, ref.cat = 1, response = "mult", ret.characteristics = F, mc.int = F){
##
## Get a design matrix from f.make.design by setting the right haplin parameters
## (Is mostly just to build a temporary sufficient info object from the parameters)
#
info <- list()
info$haplos$ref.cat <- ref.cat
info$model$xchrom <- xchrom
info$model$design <- design
info$model$poo <- poo
if(!is.null(comb.sex) && (comb.sex == "males")) info$variables$sel.sex <- 1
if(!is.null(comb.sex) && (comb.sex == "females")) info$variables$sel.sex <- 2
info$model$comb.sex <- comb.sex
info$haplos$selected.haplotypes <- rep(T, n.all)
#
# experimental
info$model$mc.int <- mc.int
info$model$hwe <- hwe
if(F){
	info$variables$covar <- 1
	info$variables$covar.codes <- 1:3
}

# Make the design matrix
.res <- f.design.make(maternal = maternal, response = response, info = info, ret.characteristics = ret.characteristics)
#
# The info object can be useful
attr(.res, "info") <- info
#
return(.res)
}
