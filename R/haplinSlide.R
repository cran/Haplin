haplinSlide <- function(filename, data, pedIndex, markers = "ALL", winlength = 1, printout = FALSE, verbose = FALSE, cpus = 1, table.output = FALSE, slaveOutfile = "", ...)
{
##
## RUN HAPLIN ON SLIDING WINDOWS
##
#.para <- "snow"
.para <- "snowfall"
#.para <- "doSMP"

## BRUKER "#i#" TIL AA KOMMENTERE VEKK DET SOM GIR ADVARSEL MED R CMD check, SIDEN DE RETTE PAKKENE IKKE ER INSTALLERT


#
## GET HAPLIN DEFAULTS, OVERRIDE DEFAULTS FOR SPECIFIED ARGUMENTS
.info <- f.catch(match.call(), formals())
#
##
if(winlength > 1) cat("\nImportant: Remember that SNPs must be in correct physical ordering\nfor a sliding window approach to be meaningful!\n")

.missdata <- missing(data)
.misspedIndex <- missing(pedIndex)

## IMPORTANT: DO NOT REMOVE MISSING YET, EVEN IF use.missing = F (I.E. REMOVE WINDOW-WISE, NOT LISTWISE)
.use.missing <- .info$model$use.missing ## SAVE OLD VALUE
.info$model$use.missing <- T ## ENFORCE






#
##


if(!.missdata && (class(data) == "gwaa.data")){
	## ONLY REDUCE NUMBER OF MARKERS, BUT DO NOT CONVERT YET
	.data.read <- data[, .info$filespecs$markers]
}else{
	## READ AND CONVERT DATA
	.data.read <- f.get.data(data, pedIndex, .info)
	.info <- attr(.data.read, "info")
}


.info$model$use.missing <- .use.missing ## REVERT TO ORIGINAL VALUE



#
## FIND MARKERS INCLUDED IN EACH WINDOW. NB: THEY NOW REFER TO MARKERS IN THE _REDUCED_ FILE, NOT THE ORIGINAL ONE
.slides <- f.windows(markers = seq(along = .info$filespecs$markers), winlength = winlength)
## USE ORIGINAL MARKER SPECIFICATION AS NAMES
.names <- .info$filespecs$markers[.slides]
.names <- matrix(.names, ncol = ncol(.slides))
.names <- f.create.tag(.names, sep = "-")
#
## PREPARE TO RUN ON EACH WINDOW
.nres <- dim(.slides)[1]
#
## REPRODUCE HAPLIN ARGUMENTS FROM .info
.args <- f.args.from.info(.info)
#
##
cat("\n")
#
## FUNCTION TO RUN HAPLIN ON A SINGLE WINDOW
.f.run <- function(i){
	cat("Running Haplin on Window '", .names[i], "' (", i, "/", .nres, ")...:  ", sep = "")
	args_ <- .args
	#
	if(cpus > 1){
		## HAPLIN MUST BE MADE AVAILABLE IN EACH RUN
		suppressPackageStartupMessages({
			require(Haplin, quietly = T)
		})
	}
	#
	##
	if(.missdata){
		if(!.info$filespecs$database){
			args_$markers <- .slides[i,]
			args_$filename <- NULL
			args_$data <- .data.read
		}else{
			args_$markers <- .info$filespecs$markers[.slides[i,]]
		}
	}else{
		args_$markers <- .slides[i,]
		if((class(data) == "gwaa.data") & !.misspedIndex){
			args_$pedIndex <- pedIndex
		}
		args_$filename <- NULL
		args_$data <- .data.read
	}
		#
		## RUN HAPLIN
		.res <- try(do.call("haplin", args_), silent = T)
		#
		## CHECK IF ERRORS
		if(class(.res) == "try-error") cat("RUN FAILED\n")
		else{
			cat("done\n")
			if(table.output) .res <- haptable(.res)
		}		
		return(.res)
}
#
## DO THE RUN
## EITHER IN SEQUENCE (SINGLE CPU) 
## OR IN PARALLEL (MULTIPLE CPUs)
if(cpus == 1){
	.res.list <- lapply(seq(length.out = .nres), .f.run)
	names(.res.list) <- .names
	#
	## CHECK FOR HAPLIN FAILURES
	.errs <- sapply(.res.list, class) == "try-error"
}else{
	#
	if(!is.numeric(cpus)) stop('The number of cpu-s "cpus" must be numeric!', call. = F)
	#
	#	## INITIALIZE CPUS
	if(.para == "snowfall"){
		## SNOWFALL
		sfInit(parallel = T, slaveOutfile = slaveOutfile, cpus = cpus)
		on.exit(sfStop())
	}
	if(.para == "doSMP"){
		## doSMP BACKEND TO foreach
#i#		w <- startWorkers(workerCount = cpus)
#i#		on.exit(stopWorkers(w))
#i#		registerDoSMP(w)
	}
	if(.para == "snow"){
		## SNOW
		# w <- makeCluster(cpus, type = type)
#i#		w <- makeCluster(cpus)
#i#		on.exit(stopCluster(w))
	}
	#
	## RUN BY SPLITTING ON CPUS
	cat("Run Haplin using ", cpus, " cpu's\n", sep = "")
	if(.para == "snowfall")	.res.list <- sfLapply(seq(length.out = .nres), .f.run)
#i#	if(.para == "snow") .res.list <- parLapply(w, seq(length.out = .nres), .f.run)
#i#	if(.para == "doSMP") .res.list <- foreach(i = seq(length.out = .nres)) %dopar% .f.run(i)
	names(.res.list) <- .names
	#
	## CHECK FOR HAPLIN FAILURES
	if(.para == "snowfall") .errs <- sfSapply(.res.list, class) == "try-error"
#i#	if(.para == "snow") .errs <- parLapply(w, .res.list, class) == "try-error"
#i#	if(.para == "doSMP") .errs <- foreach(i = seq(along = .res.list)) %dopar% class(.res.list[[i]]) == "try-error"
}
#
## GO THROUGH POSSIBLE ERRORS
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
#
## SET CLASS
if(!table.output) class(.res.list) <- "haplinSlide"
#
## RETURN RESULT LIST
return(.res.list)
}
