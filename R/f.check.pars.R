f.check.pars <- function(mcall, defaults){
##
## SET DEFAULTS AND CHECK THE LOGIC OF THE INPUT PARAMETERS TO HAPLIN
## NOTE: DOES NOT CHECK FOR SUPERFLUOUS ARGUMENTS
##
#
## REPLACE DEFAULTS WITH SPECIFIED PARAMETERS
params <- defaults
params[names(mcall)] <- mcall
#
## COLLECT IN LIST
#.info <- list(filename = info$filename, filespecs = .filespecs, model = .model, variables = .variables, haplos = .haplos, control = .control)
.info <- list(filename = 
	params[["filename"]], 
	filespecs = params[c("markers", "n.vars", "sep", "allele.sep", "na.strings", "subset", "filehash")], 
	model = params[c("design", "use.missing", "xchrom", "maternal", "test.maternal", "scoretest")], 
	variables = params[c("ccvar", "covar", "sex", "comb.sex")],
	haplos = params[c("reference", "response", "threshold", "max.haplos", "haplo.file")],
	control = params[c("resampling", "max.EM.iter", "data.out", "verbose", "printout")]
	)
class(.info) <- "info"
#
## DEFINE SPECIFIC VARIABLES
.xchrom <- .info$model$xchrom
.ccdesign <- .info$model$design %in% c("cc.triad", "cc")
#
## BASIC RECURSIVE CHECK FOR CORRECT ARGUMENT VALUES. 
## NOTE: DOES NOT CHECK NULL VALUES
.allowed <- .info
.allowed$model$scoretest <- c("yes", "no", "only")
.allowed$model$design <-  c("triad", "cc.triad", "cc")
.allowed$haplos$response <-  c("mult", "free")
.allowed$control$data.out <-  c("no", "basic", "prelim", "null", "full")
.allowed$control$resampling <- c("no", "jackknife")
.allowed$variables$comb.sex <-  c("males", "females", "single", "double")
#
for(i in seq(along = .info)[-1]){# NO CHECK FOR filename
	for(j in seq(along = .info[[i]])){
		if(length(.info[[c(i,j)]]) == 1){# EXCLUDES NULL VALUES AND VECTOR VALUES (LIKE markers)
			if(!is.element(.info[[c(i,j)]], .allowed[[c(i,j)]])){
				stop(paste("The argument ", names(.info[[i]])[j], " has an invalid value. \nIt should be one of: \n", paste(.allowed[[c(i, j)]], collapse = ", "), sep = ""), call. = F)
			}
		}
	}
}
#
##
if(.info$model$design %in% c("triad", "cc.triad")) {
	.info$model$fam <- "mfc"
}
if(.info$model$design == "cc"){
	.info$model$fam <- "c"
}
.info$filespecs$database <- F
if(!is.name(.info$filename)){
	## TEST WHETHER filename SPECIFICES A FILE OR A DIRECTORY
	.filetest.f <- file_test("-f", .info$filename)
	.filetest.d <- file_test("-d", .info$filename)
	#
	if(!.filetest.f & !.filetest.d) stop(paste(.info$filename, " is not a file nor a directory.", sep = ""), call. = F)
	if(.filetest.d){
		## ROUGH CHECK THAT IT IS A VALID DATABASE
		.info.name <- paste(.info$filename, "/000info.RData", sep = "")
		.test <- file_test("-f", .info.name)
		if(!.test) stop(paste(.info$filename, " does not look like a proper Haplin database.", sep = ""), call. = F)
		.info$filespecs$database <- T
		#
		## LOAD OLD INFO-VALUES AND USE THEM, QUIETLY IGNORE THE NEW ONES
		###.info.old <- local({ ### KOMMENTERT VEKK DETTE FOR AA UNNGAA FEILMELDING I R CMD check. BRUKER DET NEDENFOR I STEDET
		###		load(.info.name)
		###		info
		###})
		.info.old <- list()
		.info.old <- within(.info.old, load(.info.name))$info ## HAR IKKE SJEKKET DENNE....
		.info$filespecs[c("n.vars", "sep", "allele.sep")] <- .info.old$filespecs[c("n.vars", "sep", "allele.sep")]
		.info$model$design <- .info.old$model$design
	}
}
#
## SET VALUES IN THE gwaa.data CASE
if(class(mcall$data) == "gwaa.data"){
	if(.info$model$design == "triad"){
		.info$filespecs$n.vars <- 10
		if(.xchrom){
			.info$variables$sex <- 7
		}
	}
	if(.info$model$design == "cc.triad"){
		.info$filespecs$n.vars <- 10
		.info$variables$ccvar <- 10
		if(.xchrom){
			.info$variables$sex <- 7
		}
	}
	if(.info$model$design == "cc"){
		.info$filespecs$n.vars <- 4
		.info$variables$ccvar <- 4
		if(.xchrom){
			.info$variables$sex <- 3
		}
	}
}

#
## ONLY USE response = "mult" FOR MALES ONLY IN xchrom
if(.xchrom && !is.null(.info$variables$comb.sex) && .info$variables$comb.sex == "males" && !is.element(.info$haplos$response, c("mult"))){
	warning('Can only use response = "mult" with comb.sex = "males". Has been changed to "mult".', call. = F)
	.info$haplos$response <- "mult"
}
#
## TEST FOR CC VARIABLE SPECIFICATION
if(.ccdesign){
	if(is.null(.info$variables$ccvar)) stop('Parameter "ccvar" must be specified when using design "cc.triad" or "cc"!', call. = F)
	if(.info$filespecs$n.vars == 0) stop('Parameter "n.vars" must be specified when using design "cc.triad" or "cc"!', call. = F)
	if(.info$variables$ccvar > .info$filespecs$n.vars) stop('Parameter "n.vars" must be at least as large as parameter "ccvar"!', call. = F)
}else{
	if(!is.null(.info$variables$ccvar))stop('Parameter "ccvar" should only be specified when using design "cc.triad" or "cc"!', call. = F)
}
#
## MAKE SURE sex COLUMN IS SPECIFIED IF ON THE xchrom
if(!is.logical(.xchrom))stop('Argument "xchrom" must be a logical (either "TRUE" or "FALSE")', call. = F)
if(.xchrom & (.info$filespecs$n.vars == 0)) stop('Argument "n.vars" must be at least 1 to allow for a sex variable when "xchrom = TRUE"', call. = F)
if(.xchrom & !is.numeric(.info$variables$sex)) stop('Argument "sex" should be a numeric value (the column number of the sex variable) when "xchrom = TRUE"', call. = F)
if(.xchrom & (.info$model$design == "cc")) stop('The "cc" design is not (yet) implemented for x-chromosomes', call. = F)
#if(.xchrom & .info$model$maternal) stop('"xchrom = TRUE" cannot currently be used with "maternal = TRUE"', call. = F)
#if(!.xchrom & (.info$haplos$response == "xinact")) stop('Response "xinact" is only intended for the X-chromosome!', call. = F)
#
if(F){
	## CHECK THAT comb.sex IS NOT SPECIFIED UNLESS xchrom = T
	if(!.xchrom & !is.null(.info$variables$comb.sex)) warning('Argument "comb.sex" is only implemented for models where "xchrom = TRUE"', call. = F)
}
#
## SET VARIABLE sel.sex SO THAT DON'T HAVE TO REPLACE WITH comb.sex ALL THE WAY THROUGH
if(identical(.info$variables$comb.sex, "males")) .info$variables$sel.sex <- 1
if(identical(.info$variables$comb.sex, "females")) .info$variables$sel.sex <- 2
#
## MAKE SURE FULL DATA'S NOT DUMPED WHEN scoretest = "only"
if(.info$model$scoretest == "only" & (.info$control$data.out == "full")){
	warning('Since data.out = "full", scoretest argument is changed from "only" to "no"', call. = F)
	.info$control$scoretest <- "no" 
}
#
## CHECK VALUES FOR MARKERS
if(!is.numeric(.info$filespecs$markers) & !identical(.info$filespecs$markers, "ALL")) stop('"markers" argument must be either "ALL" (default) or an integer value.', call. = F) 
#
##
#
## ONLY USE response = "mult" FOR cc (FOR NOW), AND NO MATERNAL EFFECTS
if((.info$model$design == "cc") & (.info$haplos$response != "mult")){
	warning('Can only use response = "mult" with design = "cc" (for now...). Has been changed to "mult".', call. = F)
	.info$haplos$response <- "mult"
}
if((.info$model$design == "cc") & (.info$model$maternal)) stop('Cannot use maternal = TRUE with design = "cc"', call. = F)
#
## CHECK THAT response IS RESTRICTED
if((.info$haplos$response %in% c("mult")) & is.element(.info$haplos$reference, c("reciprocal", "population"))){
	warning('response = "mult" must be used with reference category (numeric or "ref.cat"). Has been changed to reference = "ref.cat"', call. = F)
	.info$haplos$reference <- "ref.cat"
}
#
## IF test.maternal IS TRUE, MAKE SURE maternal IS ALSO TRUE
if(.info$model$test.maternal) .info$model$maternal <- TRUE


## TEST OM covar, ccvar etc ER MINDRE ENN ELLER LIK n.vars, om de er heltall etc.

return(.info)


}
