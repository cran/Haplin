haplinSlide <- function(filename, 
#markers = "ALL", n.vars = 0, sep = " ", allele.sep = ";", na.strings = "NA",
#design = "triad", use.missing = FALSE, xchrom = FALSE, maternal = FALSE, test.maternal = FALSE, scoretest = "no",
#ccvar = NULL, covar = NULL, sex = NULL,
#reference = "reciprocal", response = "free", threshold = 0.01, max.haplos = NULL, haplo.file = NULL,
#resampling = FALSE, max.EM.iter = 50, data.out = FALSE, verbose = TRUE, printout = TRUE,
markers = "ALL",
winlength = 1, ...)
{
#
.mcall <- lapply(match.call()[-1], function(x) eval.parent(x, 3))
#.defaults <- formals()
.defaults <- formals(haplin)
.info <- f.check.pars(.mcall, .defaults)
#
## READ DATA
if(.info$control$verbose)	cat("\nReading data from file...  ")
if((.info$model$design == "triad") | (.info$model$design == "cc.triad")) {
	.fam <- "mfc"
}
if(.info$model$design == "cc") .fam <- "c"
.data.read <- f.read.data(indata = .info$filename, sep = .info$filespecs$sep, allele.sep = .info$filespecs$allele.sep, na.strings = .info$filespecs$na.strings, markers = .info$filespecs$markers, use.missing = .info$model$use.missing, variables = .info$filespecs$n.vars, family = .fam) ##
if(.info$control$verbose)	cat("Done\n")

.markers <- attr(.data.read, "markers")

cat("remember: SNPs must be in correct order!\n")

.slides <- f.windows(markers = .markers, winlength = winlength)
.names <- f.create.tag(.slides, sep = "-")

.args <- .mcall ## HMMM, FANGER VEL IKKE OPP EVENT. ENDRINGER I .info

.args$winlength <- NULL # TO AVOID haplin PROPER CHOKING

.nres <- dim(.slides)[1]
.res.list <- vector(.nres, mode = "list")
names(.res.list) <- .names

for(i in seq(length.out = .nres)){
	.args$markers <- .slides[i,]
	.res.list[[i]] <- do.call("haplin", .args)
}

class(.res.list) <- "haplinSlide"

return(.res.list)

}
