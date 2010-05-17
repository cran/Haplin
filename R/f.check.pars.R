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
	filespecs = params[c("markers", "n.vars", "sep", "allele.sep", "na.strings")], 
	model = params[c("design", "use.missing", "xchrom", "maternal", "test.maternal", "scoretest")], 
	variables = params[c("ccvar", "covar", "sex")],
	haplos = params[c("reference", "response", "threshold", "max.haplos", "haplo.file")],
	control = params[c("resampling", "max.EM.iter", "data.out", "verbose", "printout")]
	)
class(.info) <- "info"
#
## DEFINE SPECIFIC VARIABLES
.ccvar <- .info$variables$ccvar
.n.vars <- .info$filespecs$n.vars
.xchrom <- .info$model$xchrom
#
## BASIC RECURSIVE CHECK FOR CORRECT ARGUMENT VALUES. 
## NOTE: DOES NOT CHECK NULL VALUES
.allowed <- .info
.allowed$model$scoretest <- c("yes", "no", "only")
.allowed$model$design <-  c("triad", "cc.triad", "cc")
.allowed$haplos$response <-  c("mult", "free")
.allowed$control$data.out <-  c("no", "basic", "prelim", "null", "full")
#
for(i in seq(along = .info)[-1]){# NO CHECK FOR filename
	for(j in seq(along = .info[[i]])){
		if(length(.info[[c(i,j)]]) == 1){# EXCLUDES NULL VALUES AND VECTOR VALUES (LIKE markers)
			if(!is.element(.info[[c(i,j)]], .allowed[[c(i,j)]])){
				stop(paste("The argument ", names(.info[[i]])[j], " has an invalid value. \nIt should be one of: \n", paste(.allowed[[c(i, j)]], collapse = ", "), sep = ""))
			}
		}
	}
}
#
## TEST FOR CC VARIABLE SPECIFICATION
.ccdesign <- .info$model$design %in% c("cc.triad", "cc")
if(.ccdesign){
	if(is.null(.ccvar)) stop('Parameter "ccvar" must be specified when using design "cc.triad" or "cc"!')
	if(.n.vars == 0) stop('Parameter "n.vars" must be specified when using design "cc.triad" or "cc"!')
	if(.ccvar > .n.vars) stop('Parameter "n.vars" must be at least as large as parameter "ccvar"!')
}else{
	if(!is.null(.ccvar))stop('Parameter "ccvar" should only be specified when using design "cc.triad" or "cc"!')
}
#
## MAKE SURE sex COLUMN IS SPECIFIED IF ON THE xchrom
if(!is.logical(.xchrom))stop('Argument "xchrom" must be a logical (either "TRUE" or "FALSE")')
if(.xchrom & (.info$filespecs$n.vars == 0)) stop('Argument "n.vars" must be at least 1 to allow for a sex variable when "xchrom = TRUE"')
if(.xchrom & !is.numeric(.info$variables$sex)) stop('Argument "sex" should be a numeric value (the column number of the sex variable) when "xchrom = TRUE"')
if(.xchrom & .ccdesign) stop('The designs "cc" and "cc.triad" are not (yet) implemented for x-chromosomes')
if(.xchrom & .info$model$maternal) stop('"xchrom = TRUE" cannot currently be used with "maternal = TRUE"')
#
## MAKE SURE FULL DATA'S NOT DUMPED WHEN scoretest = "only"
if(.info$model$scoretest == "only" & (.info$control$data.out == "full")){
	warning('Since data.out = "full", scoretest argument is changed from "only" to "no"')
	.info$control$scoretest <- "no" 
}
#
## CHECK VALUES FOR MARKERS
if(!is.numeric(.info$filespecs$markers) & !identical(.info$filespecs$markers, "ALL")) stop('"markers" argument must be either "ALL" (default) or an integer value.') 
#
##
#
## ONLY USE response = "mult" FOR cc (FOR NOW), AND NO MATERNAL EFFECTS
if((.info$model$design == "cc") & (.info$haplos$response != "mult")){
	warning('Can only use response = "mult" with design = "cc" (for now...). Has been changed.')
	.info$haplos$response <- "mult"
}
if((.info$model$design == "cc") & (.info$model$maternal)) stop('Cannot use maternal = TRUE with design = "cc"')
#
## CHECK THAT response IS RESTRICTED
#if((.info$haplos$response != "free") & .info$model$maternal) stop('response != "free" not implemented for maternal effects')

if((.info$haplos$response == "mult") & is.element(.info$haplos$reference, c("reciprocal", "population"))){
	warning('response = "mult" must be used with reference category (numeric or "ref.cat"). Currently changed to reference = "ref.cat"')
	.info$haplos$reference <- "ref.cat"
}
#
## IF test.maternal IS TRUE, MAKE SURE maternal IS ALSO TRUE
if(.info$model$test.maternal) .info$model$maternal <- TRUE


## TEST OM covar, ccvar etc ER MINDRE ENN ELLER LIK n.vars, om de er heltall etc.

return(.info)


}
