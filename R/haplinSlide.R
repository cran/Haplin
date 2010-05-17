haplinSlide <- function(filename, 
#markers = "ALL", n.vars = 0, sep = " ", allele.sep = ";", na.strings = "NA",
#design = "triad", use.missing = FALSE, xchrom = FALSE, maternal = FALSE, test.maternal = FALSE, scoretest = "no",
#ccvar = NULL, covar = NULL, sex = NULL,
#reference = "reciprocal", response = "free", threshold = 0.01, max.haplos = NULL, haplo.file = NULL,
#resampling = FALSE, max.EM.iter = 50, data.out = FALSE, verbose = TRUE, printout = TRUE,
markers = "ALL",
winlength = 1, printout = FALSE, verbose = FALSE, ...)
{
#
.mcall <- lapply(match.call()[-1], function(x) eval.parent(x, 3))
#.defaults <- formals()
#
## STANDARD haplin DEFAULTS
.defaults <- formals(haplin)
#
## SPECIFIC DEFAULTS FOR haplinSlide
.defaults.slide <- formals()
#
## OVERRIDE STANDARD haplin DEFAULTS
.defaults[names(.defaults.slide)] <- .defaults.slide
#
##
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

.args <- f.args.from.info(.info)
###.args <- .mcall ## HMMM, FANGER VEL IKKE OPP EVENT. ENDRINGER I .info
###.args$winlength <- NULL # TO PROTECT haplin PROPER FROM CHOKING


.nres <- dim(.slides)[1]
.res.list <- vector(.nres, mode = "list")
names(.res.list) <- .names

cat("\n")
for(i in seq(length.out = .nres)){
	## RUN HAPLIN ON EACH WINDOW
	cat("-----------------------\n")
	cat("Running Haplin on Window '", .names[i], "'...:\n", sep = "")
	.args$markers <- .slides[i,]
	.res.list[[i]] <- try(do.call("haplin", .args), silent = T)
	if(class(.res.list[[i]]) == "try-error") cat("RUN FAILED\n")
	else cat("...done\n")
}
#
## CHECK FOR HAPLIN FAILURES
.errs <- (sapply(.res.list, class) == "try-error")
if(any(.errs)){
	## COLLECT ERROR MESSAGES, CHANGE OUTPUT TO NA WITH ERR. MESS. AS ATTRIBUTES
	.mess <- .res.list[.errs]
	.err.res <- lapply(.mess, function(x) {
		.tmpres <- NA
		attributes(x) <- NULL # REMOVE "try-error"
		attr(.tmpres, "error.message") <- x
		return(.tmpres)
	})
	.res.list[.errs] <- .err.res
}

class(.res.list) <- "haplinSlide"

return(.res.list)

}
