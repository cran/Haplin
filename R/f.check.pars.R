f.check.pars <- function(mcall, defaults){
##
## SET DEFAULTS AND CHECK THE LOGIC OF THE INPUT PARAMETERS TO HAPLIN
##
#
## SET DEFAULT VALUES
#.filespecs <- list(markers = "ALL", n.vars = 0, sep = " ", allele.sep = ";", na.strings = "NA")
#.model <- list(design = "triad", use.missing = F, xchrom = F, maternal = F)
#.variables <- list(ccvar = NULL, covar = NULL, sex = NULL)
#.haplos <- list(reference = "reciprocal", response = "free", threshold = 0.01, max.haplos = NULL, haplo.file = NULL)
#.control <- list(resampling = F, max.EM.iter = 50, data.out = F, verbose = T, printout = T)
#
## CHECK THAT ALL SPECIFIED VARIABLE NAMES ARE VALID
#if(any(!(names(info$filespecs) %in% names(.filespecs))))stop("Invalid name(s) in filespecs input")
#if(any(!(names(info$model) %in% names(.model))))stop("Invalid name(s) in model input")
#if(any(!(names(info$variables) %in% names(.variables))))stop("Invalid name(s) in variables input")
#if(any(!(names(info$haplos) %in% names(.haplos))))stop("Invalid name(s) in haplos input")
#if(any(!(names(info$control) %in% names(.control))))stop("Invalid name(s) in control input")
#
## REPLACE DEFAULTS WITH SPECIFIED PARAMETERS
params <- defaults
params[names(mcall)] <- mcall
#.filespecs[names(info$filespecs)] <- info$filespecs
#.model[names(info$model)] <- info$model
#.variables[names(info$variables)] <- info$variables
#.haplos[names(info$haplos)] <- info$haplos
#.control[names(info$control)] <- info$control
#
## COLLECT IN LIST
#.info <- list(filename = info$filename, filespecs = .filespecs, model = .model, variables = .variables, haplos = .haplos, control = .control)
.info <- list(filename = 
	params[["filename"]], 
	filespecs = params[c("markers", "n.vars", "sep", "allele.sep", "na.strings")], 
	model = params[c("design", "use.missing", "xchrom", "maternal", "scoretest")], 
	variables = params[c("ccvar", "covar", "sex")],
	haplos = params[c("reference", "response", "threshold", "max.haplos", "haplo.file")],
	control = params[c("resampling", "max.EM.iter", "data.out", "verbose", "printout")]
	)
#
## DEFINE SPECIFIC VARIABLES
.ccvar <- .info$variables$ccvar
.n.vars <- .info$filespecs$n.vars
#
## CHECK DESIGNS
.designs <- c("triad", "cc.triad", "cc")
if(!(.info$model$design %in% .designs)) stop(paste("model$design must be one of", .designs))
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
## MAKE SURE DATA'S NOT DUMPED WHEN scoretest = "only"
if(.info$model$scoretest == "only" & .info$control$data.out) warning('data.out = TRUE overrides scoretest = "only"')
#
## ONLY USE response = "mult" FOR cc (FOR NOW), AND NO MATERNAL EFFECTS
if((.info$model$design == "cc") & (.info$haplos$response != "mult")) stop('Can only use response = "mult" with design = "cc" (for now...)')
if((.info$model$design == "cc") & (.info$model$maternal)) stop('Cannot use maternal = TRUE with design = "cc"')
#
## CHECK THAT response IS RESTRICTED
#if((.info$haplos$response != "free") & .info$model$maternal) stop('response != "free" not implemented for maternal effects')
if((.info$haplos$response != "free") & is.element(.info$haplos$reference, c("reciprocal", "population"))) stop('response != "free" must be used with reference category')

## TEST OM covar, ccvar etc ER MINDRE ENN ELLER LIK n.vars, om de er heltall etc.

return(.info)


}
