hapCovar <- function(nall, n.strata = 1, cases, controls, haplo.freq, RR, RRcm, RRcf, RRstar, RR.mat, RRstar.mat, xchrom = F, sim.comb.sex = "double", BR.girls, response = "mult", ...){	
#
#
#############################################################################
## Note: This part is almost identical to hapRun
#
.sim.maternal <- FALSE
if(!missing(RR.mat) | !missing(RRstar.mat)) .sim.maternal <- TRUE
#
.sim.poo <- FALSE
if(!missing(RRcm) | !missing(RRcf)) .sim.poo <- TRUE
#
if(.sim.poo && !missing(RR)) stop("RR cannot be present at the same time as RRcm and RRcf", call. = F)
#
## The number of possible haplotypes
.nhaplo <- prod(nall)
## The number of loci
.nloci <- length(nall)
## Relative risks for controls
.RR.controls <- rep(1,.nhaplo)
#
.missing.controls <- FALSE
if(missing(controls)) .missing.controls <- TRUE
if(.missing.controls) controls <- c(mfc=0)
#
if(xchrom & sim.comb.sex %in% c("females","males")) BR.girls <- 1
#
if((n.strata==1 && length(cases)>1) | any(sapply(cases, length)>1)) stop("Each element of list cases can only have length 1", call. = F)
if(!.missing.controls && (n.strata==1 && length(controls)>1) | any(sapply(controls, length)>1)) stop("Each element of list controls can only have length 1", call. = F)
#
## Defining arguments for each stratum
.arg <- list(nall=nall, n.strata=n.strata, cases=cases, controls=controls, haplo.freq=haplo.freq, sim.maternal=.sim.maternal, sim.poo=.sim.poo, xchrom=xchrom, sim.comb.sex=sim.comb.sex, nhaplo=.nhaplo, nloci=.nloci, n.sim=10)
if(.sim.poo) .RR.arg <- list(RRcm=RRcm, RRcf=RRcf, RRstar=RRstar)
else .RR.arg <- list(RR=RR, RRstar=RRstar)
if(.sim.maternal) .RR.arg <- c(.RR.arg,list(RR.mat=RR.mat, RRstar.mat=RRstar.mat))
if(xchrom) .RR.arg <- c(.RR.arg,list(BR.girls=BR.girls))
.arg <- c(.arg, .RR.arg)
.strat.arg <- do.call(f.hapArg,args=.arg)
.f.prob.arg <- .strat.arg[, -which(colnames(.strat.arg)%in%c("cases","controls","n.sim")), drop=FALSE]
#
## Misc tests
lapply(1:n.strata, function(x){do.call(f.hapTests, args = .strat.arg[x, ])})
#
if(xchrom & sim.comb.sex == "males") message("The males are simulated assuming no contribution from fathers to sons")
#
## Choose Haplin design according to the arguments "cases" and "controls"
.design <- sapply(1:n.strata, function(x){
	if(all(.strat.arg[,"control.mat"][[x]]==0) | .missing.controls){
		if(all(.strat.arg[,"case.design"][[x]]=="c")) stop("Only case children are given. No controls are available", call. = F)
		if(all(.strat.arg[,"case.design"][[x]]=="fc") & xchrom) stop("No controls are available", call. = F)
		.design <- "triad"
	}	
	else if(all(.strat.arg[,"case.design"][[x]]=="c") & all(.strat.arg[,"control.design"][[x]]=="c")) .design <- "cc"
	else .design <- "cc.triad"
	return(.design)
})
if(length(unique(.design))==1) .design <- unique(.design)
else stop("Unable to specify haplin design due to the combination of arguments \"cases\" and \"controls\"", call.=F)
#
## Stop if cc and xchrom. Not yet implemented
if(.design == "cc" & xchrom) stop("Design \"cc\" and xchrom is not yet implemented ", call. = F)
#
## Specify arguments to be passed to haplin
.n.vars <- 0
.ccvar <- NULL
.sex = NULL
if(xchrom) .n.vars <- .n.vars + 1
if(.design!="triad") .n.vars <- .n.vars + 1
if(.design != "triad" & !xchrom) .ccvar <- .n.vars
else if(.design == "triad" & xchrom) .sex <- .n.vars	
else if(.design != "triad" & xchrom){
	.ccvar <- .n.vars - 1
	.sex <- .n.vars
}
#
## Find reference category
if(is.list(haplo.freq)) .haplo.freq <- haplo.freq[[which.max(unlist(lapply(haplo.freq,max)))]]
else .haplo.freq <- haplo.freq
.ref.cat <- which.max(.haplo.freq)
.lu <- list(...)
if("reference" %in% names(.lu)) .ref.cat <- .lu$reference
#
## Response
.response <- response
#
## Arguments to Haplin
.haparg <- list(n.vars = .n.vars, design = .design, ccvar = .ccvar, xchrom = xchrom, sex = .sex, verbose = FALSE, use.missing = TRUE, threshold = 0, reference = .ref.cat, data.out = "prelim")
#
## Arguments for each strata
.f.prob.arg <- .f.prob.arg[, -which(colnames(.f.prob.arg)%in%c("gen.missing.cases", "gen.missing.controls", "nloci", "nhaplo", "case.des", "control.des", "case.mat", "control.mat", "case.design", "control.design", "nall", "sim.poo")), drop=FALSE]
colnames(.f.prob.arg)[which(colnames(.f.prob.arg)=="xchrom")] <- "sim.xchrom" 
if(!.sim.poo){
	.f.prob.arg.RRcmcf <- cbind(.f.prob.arg[,"RR"],.f.prob.arg[,"RR"])
	colnames(.f.prob.arg.RRcmcf) <- c("RRcm","RRcf")
	.f.prob.arg <- cbind(.f.prob.arg, .f.prob.arg.RRcmcf)
	.f.prob.arg <- .f.prob.arg[,-which(colnames(.f.prob.arg)=="RR"), drop=FALSE]
}
#################################################################
#
## Calculate the asymptotic covariance matrix for each strata
.var.covar <- list()
for(i in 1:n.strata){
	.tmp.strat.arg <- .strat.arg[i,]
	####################################
	## Beta values
	## Compute beta values for the specified relative risks
	.RR <- as.data.frame(.tmp.strat.arg[which(grepl("RR",names(.tmp.strat.arg)))])
	.RRstar <- .RR[,which(grepl("star",names(.RR))), drop=F]
	#
	if(.response=="mult"){
		if(any(apply(.RRstar,2,sum)!=.nhaplo)) stop("Arguments RRstar and/or RRstar.mat do not correspond to a multiplicative dose-response model", call. = F)
		.RR <- .RR[,-which(grepl("star",names(.RR))), drop=F]
	}
	.names <- names(.RR)
	.RR.ref <- .RR[.ref.cat,,drop=F]
	.RR <- lapply(1:ncol(.RR), function(x) .RR[,x] <- .RR[,x]/.RR.ref[,x])
	#
	.RR <- lapply(1:length(.RR), function(x){
		if(!grepl("star",.names[x])) .RR <- .RR[[x]][-.ref.cat]
		else if(grepl("star",.names[x]) & .nhaplo<=2) .RR <- .RR[[x]][-.ref.cat]
		else .RR <- .RR[[x]]
		.RR
	})
	.RR.beta <- as.vector(log(unlist(.RR)))
	#
	## Compute beta values for haplo freq
	.haplo.coef <- .tmp.strat.arg$"haplo.freq"
	.haplo.beta <- f.beta.haplo.freq.asymp(haplo.freq = .haplo.coef)
	#
	## Compute beta values for cc and xchrom
	.beta <- c(.haplo.beta,.RR.beta)
	.k <- 0
	if(.design!="triad"){
		.k <- .tmp.strat.arg$"controls"/.tmp.strat.arg$"cases"
		.beta <- c(0,.beta)
	}
	if(xchrom) .beta <- c(-log(1/BR.girls), .beta)
	
	##########################################
	## Probabilities
	
	.design.matrix <- f.design.get(n.all = .nhaplo, design = .design, xchrom = xchrom, maternal = .sim.maternal, poo = .sim.poo, hwe = T, comb.sex = sim.comb.sex, ref.cat = .ref.cat, response = .response, ret.characteristics = F, mc.int = F)
	.info <- attr(.design.matrix, "info")
	.X <- as.matrix(.design.matrix)
	names(.beta) <- colnames(.X)
	.prob <- f.prob.asymp(beta = .beta, design = .design, X = .X, k = .k)
	
	############################################
	## Grid	and data file
	## Make grid 
	.design.grid <- f.design.get(n.all = .nhaplo, design = .design, xchrom = xchrom, maternal = .sim.maternal, poo = .sim.poo, hwe = T, comb.sex = sim.comb.sex, ref.cat = .ref.cat, response = .response, ret.characteristics = T, mc.int = F)
	#
	## cc and xchrom not yet implemented
	#if((.design == "cc") & xchrom) .grid <- .design.grid
	#
	.grid <- expand.grid(lapply(.design.grid, function(x){1:x}))
	#
	.case.design = .tmp.strat.arg$"case.design"
	.control.design = .tmp.strat.arg$"control.design"
	.grid <- f.grid.asymp(pos = nrow(.grid), design = .design, xchrom = xchrom, n.vars = .n.vars, nall = nall, case.design = .case.design, control.design = .control.design)
	#
	.ncells <- nrow(.grid)
	.haparg$data <- .grid
	#
	## Compute data file prepared for analysis
	.data <- do.call("haplin0", args=.haparg)
	if(.design!="triad"){
		.data$cc[.data$cc=="case"] <- 2
		.data$cc[.data$cc=="control"] <- 1
		mode(.data$cc) <- "numeric"
	}	
	if(xchrom){
		.data$sex[.data$sex=="girl"] <- 2
		.data$sex[.data$sex=="boy"] <- 1
		mode(.data$sex) <- "numeric"
	}
	#
	.orig <- sort(unique(.data$orig.lines)) # NOTE! assumes all data files derived the same way
	if(!identical(.orig, 1:.ncells)) stop() # line numbers in grid. But this may change later if grid is different than data
	.norig <- length(.orig)
	
	#######################################################
	## Compute covariance matrix
	.var.covar.strat <- f.var.covar.asymp(X = .X, data = .data, pred = .prob, ncells = .ncells, norig = .norig, orig = .orig, info = .info)
	#
	if(xchrom | .design!="triad") .beta <- .beta[-which(names(.beta)%in%c("cc","sex"))]
	#
	## Normalize
	.var.covar.strat <- .var.covar.strat[names(.beta),names(.beta)]/(.strat.arg[i, "cases"]$cases + .strat.arg[i, "controls"]$controls)
	if(.sim.poo) .var.covar.strat <- f.post.poo.diff(list(as.matrix(.beta)),list(as.matrix(.var.covar.strat)))
	else .var.covar.strat <- f.post.diff(list(.beta),list(.var.covar.strat))
	#
	.var.covar[[i]] <- .var.covar.strat

}
#
## Combine results from all strata
.coef <- sapply(.var.covar, function(x) x$coeff)
.covar <- sapply(.var.covar, function(x) x$covar)
.asymp <- list(coef = .coef, cov = .covar)
attr(.asymp,"ref.cat") <- .ref.cat
#
## Return beta values and covariance matrix. Reference category as attribute
return(.asymp)
}