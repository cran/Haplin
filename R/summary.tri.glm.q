"summary.tri.glm"<-
function(object, reference.method, conf.int = T, n.sim = 10000, level = 0.95, design, info, ...)
{
# CREATES A SUMMARY OF AN OBJECT OF CLASS tri.glm
# NOTE: ... IS JUST IGNORED
#
#### PREPARE: ###################
	library(MASS) 
	#cat("dette....\n")
	if(missing(info)) info <- NULL #
	.n.all <- object$nall
	.n.tri <- object$ntri
	.maternal <- object$maternal
	.res <- object$result
	.ref.cat <- object$ref.cat
	.summary <- summary.glm(.res, correlation = F)
	.coef <- .summary$coefficients[, 1] # WORKS FOR BOTH SPLUS AND R
	.cov <- .summary$cov.unscaled	##
#
	.cov.correct <- attr(object, "cov.correct")
	if(!is.null(.cov.correct)) .cov <- .cov.correct
#
	.cov.resamp <- attr(object, "cov.resamp")
	if(!is.null(.cov.resamp)) .cov <- .cov.resamp
#



#### KAN VEL SIKKERT FJERNE DENNE SENERE:

.test <- coef.tri.glm(object)
if(!all.equal(.test$coef, .coef) | !all.equal(.test$cov, .cov)) warning("No problem, really...")


###
###
###	cat("PASS PAA HVILKEN MATRISE SOM BLIR VALGT!...\n")
###
###





#	
#	
#### CONFIDENCE INTERVALS NOT REQUESTED:  #######
	if(!conf.int) {
	if(design != "triad") stop("Not implemented")
#
#   COMPUTE REPARAMETRIZATION:
		.effects <- t(f.compute.effects(.res$coefficients, design = design, maternal = .maternal, reference.method = reference.method, ref.cat = .ref.cat, n.all = .n.all, info = info))
		.ut <- list(effects = .effects, design = design, pvalues = NULL, maternal = .maternal, reference.method = reference.method, conf.int = conf.int, n.all = .n.all, n.tri = .n.tri, level = level, orig.call = object$orig.call, date = object$date, ref.cat = .ref.cat)
	}
#
#

if(F){
	### SE PÅ MULIGHETEN FOR Å UTELUKKE PROBLEMESTIMATER:
	###.prob.par <- sqrt(diag(.cov)) > 1e2
	.prob.par <- sqrt(diag(.cov))/abs(.coef) > 5
	if(all(.prob.par)) stop("Too much uncertainty in estimates")
	.cov.red <- .cov
	.cov.red[.prob.par,] <- 0
	.cov.red[,.prob.par] <- 0

	#tull <<- list(mu = .coef, cov = .cov, cov.red = .cov.red)
	.cov <- .cov.red
}



#### CONFIDENCE INTERVALS REQUESTED: ############
	if(conf.int) {
#
#  SIMULATE MULTIVARIATE DATA FOR COMPUTING CONFIDENCE INTERVALS:
		f.vis("Merk: Gjor simuleringer, bruker set.seed!", vis = F)
		set.seed(24)
		.sim <- mvrnorm(n.sim, mu = .coef, Sigma = .cov)
		dimnames(.sim) <- list(NULL, names(.coef))#

#         COMPUTE REPARAMETRIZATION:
		.effects <- f.compute.effects(.sim, design = design, maternal = .maternal, reference.method = reference.method, ref.cat = .ref.cat, n.all = .n.all, info = info)
####		if(design == "triad") .effects <- f.compute.effects(.sim, maternal = .maternal, reference.method = reference.method, ref.cat = .ref.cat, n.all = .n.all)
####		if(design == "cc") .effects <- f.compute.effects.cctemp(.sim, maternal = .maternal, reference.method = reference.method, ref.cat = .ref.cat, n.all = .n.all)
		
#
		.f.quant <- function(x){
			if(any(is.na(x))){ 
				rep(as.numeric(NA), 3)
			}else{
				quantile(x, probs = c(0.5, (1 - level)/2, 1 - (1 - level)/2))
				# SENDS ANY PROBLEMS, LIKE NA, STRAIGHT THROUGH, BUT WITH WARNING
			}
		}
		if(any(is.na(.effects))) warning("NAs in confidence intervals")
		.effects.CI <- t(apply(.effects, 2, .f.quant))
		dimnames(.effects.CI)[[2]] <- c("est.", "lower", "upper")
		.pvalues <- t(apply(.effects, 2, function(x){
		.sum <- ifelse(median(x) > 1, sum(x <= 1), sum(x >= 1))
		min(.sum/length(x) * 2, 1)
		}))
#
		.effects.CI <- cbind(.effects.CI, p.value = as.numeric(.pvalues))
		.ut <- list(effects = .effects.CI, design = design, pvalues = .pvalues, maternal = .maternal, reference.method = reference.method, conf.int = conf.int, n.all = .n.all, n.tri = .n.tri, level = level, orig.call = object$orig.call, date = object$date, ref.cat = .ref.cat)
	}
		
	class(.ut) <- "summary.tri.glm"

	return(.ut)
}
