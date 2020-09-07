snpPower <- function(cases, controls, RR, MAF, alpha = 0.05){
##
## Power calculation for a single SNP. Allows for power calculations of several combinations simultaneously.
##
## cases is a list of the number of case families. Each element contains the number of families of a specified family design. The names of the elements, i.e. the family designs, must be "mfc" (full triad), "mc" (mother-child dyad), "fc" (father-child dyad) or "c" (a single case child). 
## controls is a list of the number of control families. Each element contains the number of families of a specified family design. The possible family designs are "mfc" (full triad), "mc" (mother-child-dyad), "fc" (father-child dyad), "mf" (mother-father dyad), "c" (a single control child), "m" (a single control mother) or "f" (a single control father), i.e., the name of each element must equal one of the options. 
## RR is a numeric vector of the relative risks.
## MAF is a numeric vector of the minor allele frequencies.
## alpha is a numeric vector of the Type I Errors. Equals 0.05 by default.
##
## All arguments must be of equal length or have length equal to one.
##
#
## Power function based on normal approximation on log(OR)
or.power <- function(p2, OR, n1, n2, alpha = 0.05){
	p1 <- p2*OR/(1-p2+p2*OR)  # p1 er allelefrekvens i case-gruppen (MAF i case-gruppen). p2 er MAF i kontrollgruppen
	.sd <- sqrt((1/(n1*p1*(1-p1)))+1/(n2*p2*(1-p2)))
	return(1+pnorm(-qnorm(1-alpha/2)-log(OR)/.sd)-pnorm(qnorm(1-alpha/2)-log(OR)/.sd))
}
#
## Possible family designs
.case.des <- c("mfc", "mc", "fc", "c")
.control.des <- c("mfc", "mc", "fc", "mf", "c", "m", "f")
#
## Checking if all arguments have the same length or have length equal to one
.check.length <- unique(c(sapply(cases, length), sapply(controls, length), length(RR), length(MAF), length(alpha)))
if(length(.check.length[.check.length != 1]) >= 2) stop("All arguments must be of equal length or have length equal to 1", call. = F)
#
#.l <- max(.check.length)
#
.case.mat <- do.call("cbind", cases)
.case.name <- colnames(.case.mat)
.control.mat <- do.call("cbind", controls)
.control.name <- colnames(.control.mat)

#
## Misc errors
if(!all(names(cases) %in% .case.des) || any(duplicated(names(cases)))) stop("Argument \"cases\" is not specified correctly", call. = F)
if(!is.numeric(.case.mat) || any(.case.mat < 0) || any(rowSums(.case.mat) == 0)) stop("Argument \"cases\" must contain non-negative integers, and the values cannot all equal zero.", call. = F)
if(!all(names(controls) %in% .control.des) || any(duplicated(names(controls)))) stop("Argument \"controls\" is not specified correctly", call. = F)
if(!is.numeric(.control.mat) || any(.control.mat < 0)) stop("Argument \"controls\" must contain only non-negative integers", call. = F)
if(!is.numeric(RR) || any(RR  <= 0)) stop("The relative risk must be numeric and larger than 0", call. = F)
if(!is.numeric(MAF) || any(MAF <= 0 | MAF >= 1)) stop("\"MAF\" contains invalid value(s)", call. = F)
if(!is.numeric(alpha) || any(alpha < 0 | alpha > 1)) stop("\"alpha\" contains invalid Type I Error(s)", call. = F)
#
## Weights for the different family designs
.tab.w <- data.frame(.case.mat, .control.mat, cbind(RR = RR, MAF = MAF, alpha = alpha))
#
.rat.loss.pseudo <- apply(.tab.w, 1, function(x){
	.k <- (1-x["MAF"])/(x["MAF"]*x["RR"])
	.rat.loss <- .k/(1 + .k)^2
	.rat.loss
})
#
.rat.loss.controls <- apply(.tab.w, 1, function(x){
	.k <- (1-x["MAF"])/x["MAF"]
	.rat.loss <- .k/(1 + .k)^2
	.rat.loss
})
#
.l <- length(.rat.loss.pseudo)
.case.w <- c(2, 2, 2, 2)
.pseudo.w <- matrix(c(rep(2,.l), (1-.rat.loss.pseudo), (1-.rat.loss.pseudo), rep(0,.l)), nrow=.l)
.control.w <- matrix(c(rep(4,.l), (3-.rat.loss.controls), (3-.rat.loss.controls), rep(4, .l), rep(2,.l), rep(2,.l), rep(2,.l)),nrow=.l)
#
.case.mat0 <- matrix(rep(0, .l*4), ncol = 4)
if(nrow(.case.mat)==1){
	.case.mat <- matrix(rep(.case.mat, each=.l), ncol=length(.case.mat))
	colnames(.case.mat) <- .case.name
}
colnames(.case.mat0) <- .case.des
#
.control.mat0 <- matrix(rep(0, .l*7), ncol = 7)
if(nrow(.control.mat)==1){
	.control.mat <- matrix(rep(.control.mat, each=.l), ncol=length(.control.mat))
	colnames(.control.mat) <- .control.name
}
colnames(.control.mat0) <- .control.des
#
.case.mat0[, colnames(.case.mat)] <- .case.mat
.control.mat0[, colnames(.control.mat)] <- .control.mat
#
## Number of case- and control alleles
.n1 <- rowSums(t(t(.case.mat0)*.case.w)) ## Number of case alleles
.n0.pseudo <- rowSums(.case.mat0*.pseudo.w) ## Number of control alleles from case families
.n0 <- rowSums(.control.mat0*.control.w) ## Number of control alleles from control families
.n <- data.frame(n1 = .n1, n0.pseudo = .n0.pseudo, n0 = .n0)
.n$n0 <- .n$n0.pseudo + .n$n0 ## Number of total control alleles
.n$n0.pseudo <- NULL
#
if(any(.n$n0 == 0)) stop("The power is computed for combination(s) without control alleles", call. = F)
#
## Create data frame from the desired combinations
names(cases) <- paste("cases", names(cases), sep = ".")
names(controls) <- paste("controls", names(controls), sep = ".")
.tab <- cbind(cases, controls, .n, RR = RR, MAF = MAF, alpha = alpha)
#
## Power calculations using the normal approximation on log(OR)
.or.power <- apply(.tab, 1, function(x){
	or.power(p2 = x["MAF"], OR = x["RR"], n1 = x["n1"], n2 = x["n0"], alpha = x["alpha"])
})
#
## Return output
.output <- cbind(.tab, power = .or.power)
.output[, c("n0","n1")] <- list(NULL, NULL)
return(.output)
}
