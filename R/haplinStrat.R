haplinStrat <- function(filename, covar, ...){
#
.mcall <- lapply(match.call()[-1], function(x) eval.parent(x, 3))
#.defaults <- formals()
.defaults <- formals(haplin)
.info <- f.check.pars(.mcall, .defaults)



cat("\n## Running haplinStrat ##\n")

if(F){
	## ARGUMENTS NOT GIVEN ARE SET TO NULL (THEN LATER TO DEFAULT VALUES BY f.check.pars)
	if(missing(filename)) filename <- NULL
	if(missing(filespecs)) filespecs <- NULL
	if(missing(model)) model <- NULL
	if(missing(variables)) variables <- NULL
	if(missing(haplos)) haplos <- NULL
	if(missing(control)) control <- NULL
	#
	## COLLECT INFO, SET DEFAULT VALUES AND CHECK FOR SANITY
	.info <- list(filename = filename, filespecs = filespecs, model = model, variables = variables, haplos = haplos, control = control)
	.info <- f.check.pars(.info)
}
.info0 <- .info # INITIAL VALUES, AS SPEC. BY USER AND DEFAULTS
#
design <- .info$model$design

if(F){
	#
	## MAKE SURE CONTRASTS ARE RIGHT (ALTHOUGH NOT STRICTLY NECESSARY):
	.old.options <- options()
	on.exit(options(.old.options))
	options(contrasts = c("contr.treatment", "contr.poly"))	#
	#
	## INSTALL MASS (FOR THE mvrnorm FUNCTION):
	require(MASS)
	#	
}

#
## READ DATA AND REPORT MISSING:
cat("\nReading data from file...")
if((design == "triad") | (design == "cc.triad")) {
	.fam <- "mfc"
}
if(design == "cc") .fam <- "c"
.data.read <- f.read.data(indata = .info$filename, sep = .info$filespecs$sep, allele.sep = .info$filespecs$allele.sep, na.strings = .info$filespecs$na.strings, markers = "ALL", use.missing = .info$model$use.missing, variables = .info$filespecs$n.vars, family = .fam) ##
cat("Done\n")


.data <- .data.read

#print(head(.data))

#
## AD HOC STUFF, TO BE CHANGED:
#.data.var <- .data[, seq(length.out = .info$filespecs$n.vars), drop = F]
#.data.gen <- .data[, -seq(length.out = .info$filespecs$n.vars), drop = F]
#.ind <- matrix(1:(dim(.data.gen)[2]), nrow = 2)
#.tmp <- vector(dim(.ind)[2], mode = "list")
#for(i in seq(length.out = dim(.ind)[2])){
#    .tmp[[i]] <- paste(.data.gen[,.ind[1,i]], .data.gen[,.ind[2,i]], sep = ";")
#}
#.data.gen <- do.call("data.frame", .tmp)
#names(.data.gen) <- NULL
#.data <- data.frame(.data.var, .data.gen)

#
### DEFINE STRAT VARIABLE AND DISPLAY FREQUENCY DISTRIBUTION
.strata <- .data.read[, .info$variables$covar]
cat("\nFrequency distribution of selected stratification variable:\n")
.tmp <- table(.strata)
names(dimnames(.tmp)) <- NULL
print(.tmp)
.strata.list <- sort(unique(.strata))

.ut.list <- vector(length(.strata.list) + 1, mode = "list")
names(.ut.list) <- c("all", .strata.list)


.tmpfilename <- tempfile(tmpdir = ".")
.tmphaplofile <- tempfile(tmpdir = ".")
on.exit({unlink(.tmpfilename); unlink(.tmphaplofile)})

write.table(.data, file = .tmpfilename, sep = "\t", col.names = F, row.names = F, quote = F)

cat("\nRunning Haplin for full data file...")
###.parstmp <- .info0
.parstmp <- .mcall ## HMMM, FANGER VEL IKKE OPP EVENT. ENDRINGER I .info
if(F){
	.parstmp$filespecs$sep <- "\t"
	.parstmp$filespecs$allele.sep <- "\t"
	.parstmp$filename <- .tmpfilename
	.parstmp$control$verbose <- F
	.parstmp$control$printout <- F
}
.parstmp$sep <- "\t"
.parstmp$allele.sep <- "\t"
.parstmp$filename <- .tmpfilename
.parstmp$verbose <- F
.parstmp$printout <- F

.ut.list[["all"]] <- do.call("haplin", .parstmp)
cat("Done\n")


write.table(haptable(.ut.list[["all"]]), file = .tmphaplofile, quote = F)
#tull <<- .tmphaplofile

#return(.ut.list[["all"]])

#stop()
for(i in seq(along = .strata.list)){
	cat('\nRunning Haplin for stratum "', .strata.list[i], '"...', sep = "")
	.tmpd <- .data[.strata == .strata.list[i], ]
	write.table(.tmpd, file = .tmpfilename, sep = "\t", col.names = F, row.names = F, quote = F)
	###.parstmp <- .info0
	.parstmp <- .mcall
	if(F){
		.parstmp$filename <- .tmpfilename
		.parstmp$control$verbose <- F
		.parstmp$control$printout <- F
		.parstmp$haplos$haplo.file <- .tmphaplofile
		.parstmp$haplos$reference <- .ut.list[["all"]]$info$haplos$ref.cat
		.parstmp$filespecs$allele.sep <- "\t"
	#	.parstmp$filespecs$markers <- "ALL"
	#	.ut.list[[as.character(.strata.list[i])]] <- try(do.call("haplin", .parstmp)) # FUNKER, MEN kan RISIKERE AT NOEN BRUKER NAVNET "all" I FAKTOR
		}
		.parstmp$filename <- .tmpfilename
		.parstmp$verbose <- F
		.parstmp$printout <- F
		.parstmp$haplo.file <- .tmphaplofile
		.parstmp$reference <- .ut.list[["all"]]$info$haplos$ref.cat
		.parstmp$allele.sep <- "\t"

	.ut.list[[i + 1]] <- try(do.call("haplin", .parstmp))
#	.ut.list[[i]] <- do.call("haplin", .parstmp)
	cat("Done\n")
}

class(.ut.list) <- "haplinStrat"
return(.ut.list)


}
