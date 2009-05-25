"f.compute.effects"<-
function(x, maternal, reference.method, ref.cat, n.all, design, info)
{
# COMPUTES ALLELE FREQUENCIES, RELATIVE RISKS ETC. USING THE RESULTS
# FROM f.tri.glm. x IS THE VECTOR OF ESTIMATED COEFFICIENTS, OR POSSIBLY
# A MATRIX WHERE EACH ROW IS A SET OF ESTIMATED COEFFICIENTS
#
#### PREPARE: ##########################
.safe <- T # USES A MORE STABLE (BUT MORE TIME-CONSUMING) COMPUTATION.
		# SHOULD NOT BE NECESSARY IF NONE OF THE ALLELES ARE RARE
if(!is.matrix(x)) {
	.resmat <- matrix(x, nrow = 1)
	dimnames(.resmat) <- list(NULL, names(x))
}
else .resmat <- x #
#
#cat("dette....\n")
#if(missing(info)) .resp <- NULL
#else .resp <- info$haplos$response
#
if(design == "cc" | design == "cc.triad") .resmat <- .resmat[,-1] ## IN CASE-CONTROL DESIGN, REMOVE COLUMN FOR CASE LEVEL
#
## FILL IN ZEROS FOR NON-ESTIMATED EFFECTS
.resmat <- f.fill.effects(resmat = .resmat, info = info)
#
if(F){
.test <- f.fill.effects(resmat = .resmat, info = info)
## SHOULD NO LONGER BE NEEDED
if(n.all == 2){
	#.mf <- .resmat[,1:2]
	#.c <- .resmat[,3, drop = F]
	#.cdd <- .resmat[,4, drop = F]
	#if(maternal){
	#    .m <- .resmat[,5, drop = F]
	#    .mdd <- .resmat[,6, drop = F]
	#}else{
	#    .m <- NULL
	#    .mdd <- NULL
	#}
	
	if(!maternal){
		if(ref.cat == 1) {
			.resmat <- cbind(.resmat[, 1:2, drop = F], 0, .resmat[, 3, drop = F], 0, .resmat[, 4, drop = F])
			dimnames(.resmat)[[2]][c(3, 5)] <- paste(c("c", "cdd"), ref.cat, sep = "")
		}
	#
		if(ref.cat == 2) {
			.resmat <- cbind(.resmat[, 1:3, drop = F], 0, .resmat[, 4, drop = F], 0)
			dimnames(.resmat)[[2]][c(4, 6)] <- paste(c("c", "cdd"), ref.cat, sep = "")
		}
	}# END if(!maternal)
	#
	if(maternal){
		if(ref.cat == 1) {
			.resmat <- cbind(.resmat[, 1:2, drop = F], 0, .resmat[, 3, drop = F], 0, .resmat[, 4, drop = F], 0, .resmat[, 5, drop = F], 0, .resmat[, 6, drop = F])
			dimnames(.resmat)[[2]][c(3, 5, 7, 9)] <- paste(c("c", "cdd", "m", "mdd"), ref.cat, sep = "")
		}
	#
		if(ref.cat == 2) {
			.resmat <- cbind(.resmat[, 1:3, drop = F], 0, .resmat[, 4, drop = F], 0, .resmat[, 5, drop = F], 0, .resmat[, 6, drop = F], 0)
			dimnames(.resmat)[[2]][c(4, 6, 8, 10)] <- paste(c("c", "cdd", "m", "mdd"), ref.cat, sep = "")
		}
	}# END if(maternal)
}# END if(n.all == 2)

	if(n.all >= 3 & !maternal){# INSERT ZERO-COLUMN FOR THE REFERENCE CATEGORY (SOMEWHAT ROUNDABOUT...)
		.ind <- seq(length = dim(.resmat)[[2]])
		.insert.pos <- n.all + ref.cat
		.ind <- f.insert.vector(target = .ind, insert = NA, pos.target = .insert.pos, len.insert = 1)
		.resmat <- .resmat[,.ind, drop = F]
		dimnames(.resmat)[[2]][.insert.pos] <- paste("c", ref.cat, sep = "")
		.resmat[,.insert.pos] <- 0
		f.vis(.resmat[1,, drop = F], vis = F)
		
	}	#
#
	if(n.all >= 3 & maternal){# INSERT ZERO-COLUMNS FOR THE REFERENCE CATEGORY
		.ind <- seq(length = dim(.resmat)[[2]])
		.insert.pos <- c(n.all + ref.cat, 3 * n.all - 1 + ref.cat)
		.ind <- f.insert.vector(target = .ind,insert = rep(NA, 2), pos.target = .insert.pos, len.insert = c(1,1))
		.resmat <- .resmat[,.ind, drop = F]
		.insert.pos <- .insert.pos + c(0,1)
		dimnames(.resmat)[[2]][.insert.pos] <- paste(c("c", "m"), ref.cat, sep = "")
		.resmat[,.insert.pos] <- 0
		f.vis(.resmat[1,, drop = F], vis = F)
		
	}	#
#
	f.vis(head(.resmat), vis = T)
	f.vis(head(.test), vis = T)
	if(!identical(.resmat, .test))stop("ooops")
}## END SHOULD NO LONGER BE NEEDED

#
.J <- matrix(1, nrow = n.all, ncol = 1) #
#
#### COMPUTE ALLELE FREQUENCIES: #############
.p <- exp(.resmat[, 1:n.all, drop = F])
.p <- .p/as.numeric(.p %*% .J)	#
#
#### COMPUTE RELATIVE RISKS FOR CHILD: #############
	.R <- exp(.resmat[, (n.all + 1):(2 * n.all), drop = F])	#
	.Rstar <- exp(.resmat[, (2 * n.all + 1):(3 * n.all), drop = F])	#
# IN REFCAT PARAMETRIZATION:
		.Rstartilde <- .R^2 * .Rstar
# IN POPULATION REFERENCE:
	if(reference.method == "population") {
		.B <- as.numeric((.R * .p) %*% .J)	# POPULATION BASELINE
		.Rtilde <- .R/.B
		.Rstartilde <- .Rstartilde/.B
	}
# IN RECIPROCAL REFERENCE:
	if(reference.method == "reciprocal"){ #
#
		.f.minus <- function(x) as.numeric(x %*% rep(1, dim(x)[2])) - x #
#
		.f.kryssprod <- function(x, f.minus){
			(x * f.minus(x)) %*% rep(1, dim(x)[2])		
		} #
#
		.f.kryssprod.minus <- function(x)
			{
			.d2 <- dim(x)[2]
			(as.numeric(x %*% rep(1, .d2)) - x)^2 - (as.numeric(x^2 %*% rep(1, .d2)) - x^2)
			} #
#
		.f.minus.safe <- function(x){
			.d <- dim(x)
			.ut <- x
			.ut[,] <- NA
			for(i in seq(length = .d[2]))
			{
				.ut[,i] <- x[, -i, drop = F] %*% rep(1, .d[2] -1)
			}
			.ut
		} #
#		
		.f.kryssprod.minus.safe <- function(x, f.minus, f.kryssprod){
			.d <- dim(x)
			.ut <- x
			.ut[,] <- NA
			for(i in seq(length = .d[2]))
			{
			.ut[,i] <- f.kryssprod(x[,-i, drop = F], f.minus)
		}
			.ut
		} #
#
#
		.pR <- .p * .R #
#
		if(.safe){
		.now <- proc.time()[1]
		.Pi.safe <- .R * .f.minus.safe(.pR) / .f.minus.safe(.p)
		.Pi.minus.safe <- .f.kryssprod.minus.safe(.pR, .f.minus.safe, .f.kryssprod) / .f.kryssprod.minus.safe(.p, .f.minus.safe, .f.kryssprod)
		.Fi.safe <- .Pi.safe/.Pi.minus.safe
		.Fi.startilde.safe <- .Rstartilde/.Pi.minus.safe
		f.vis(paste("used", proc.time()[1] - .now, "\n"), vis = F)
		.Fi <- .Fi.safe
		.Fi.startilde <- .Fi.startilde.safe
}# END if(.safe)
		if(!.safe){
		.now <- proc.time()[1]
		.Pi <- .R * .f.minus(.pR) / .f.minus(.p)
##		.Pi.minus <- (.f.minus(.pR)^2 - .f.minus(.pR^2))/(.f.minus(.p)^2 - .f.minus(.p^2))
		.Pi.minus <- .f.kryssprod.minus(.pR)/.f.kryssprod.minus(.p)
		if(any(.Pi.minus == 0)){
			.Pi.minus[.Pi.minus == 0] <- 1e-10  
			warning("NAs generated during simulation, perhaps due to rare haplotype!")
		}
		.Fi <- .Pi/.Pi.minus
		.Fi.startilde <- .Rstartilde/.Pi.minus
		f.vis(paste("used", proc.time()[1] - .now, "\n"), vis = F)
}# END if(!.safe)
	}# END reciprocal
#
#
if(maternal) {#
#### COMPUTE RELATIVE RISKS FOR MATERNAL EFFECTS, IF APPLICABLE: #############
	.Rm <- exp(.resmat[, (3 * n.all + 1):(4 * n.all), drop = F])	#
	.Rmstar <- exp(.resmat[, (4 * n.all + 1):(5 * n.all), drop = F]) # 
# IN REFCAT PARAMETRIZATION:
		.Rmstartilde <- .Rm^2 * .Rmstar
# IN POPULATION REFERENCE:
	if(reference.method == "population") {
		.Bm <- as.numeric((.Rm * .p) %*% .J)	# POPULATION BASELINE
		.Rmtilde <- .Rm/.Bm
		.Rmstartilde <- .Rmstartilde/.Bm
	}
# IN RECIPROCAL REFERENCE:
	if(reference.method == "reciprocal"){
		.pRm <- .p * .Rm #
#
		if(.safe){
		.now <- proc.time()[1]
		.Pim.safe <- .Rm * .f.minus.safe(.pRm) / .f.minus.safe(.p)
		.Pim.minus.safe <- .f.kryssprod.minus.safe(.pRm, .f.minus.safe, .f.kryssprod) / .f.kryssprod.minus.safe(.p, .f.minus.safe, .f.kryssprod)
		.Fim.safe <- .Pim.safe/.Pim.minus.safe
		.Fim.startilde.safe <- .Rmstartilde/.Pim.minus.safe
		f.vis(paste("used", proc.time()[1] - .now, "\n"), vis = F)
		.Fim <- .Fim.safe
		.Fim.startilde <- .Fim.startilde.safe
}# END SAFE
		if(!.safe){
		.now <- proc.time()[1]
		.Pim <- .Rm * .f.minus(.pRm) / .f.minus(.p)
		.Pim.minus <- .f.kryssprod.minus(.pRm)/.f.kryssprod.minus(.p)
		if(any(.Pim.minus == 0)){
			.Pim.minus[.Pim.minus == 0] <- 1e-10  
			warning("NAs generated during simulation, perhaps due to rare haplotype!")
		}
		.Fim <- .Pim/.Pim.minus
		.Fim.startilde <- .Rmstartilde/.Pim.minus
		f.vis(paste("used", proc.time()[1] - .now, "\n"), vis = F)
}# END if(!.safe)
}# END RECIPROCAL
}# END if(maternal)
#
#
#
#### OUTPUT: ##########################  
	if(reference.method == "ref.cat" & !maternal) {
		#if(design == "cc"){
		#.ut <- cbind(.p, .R, .Rstar)
		#cat("gl\n")
		#} else
		.ut <- cbind(.p, .R, .Rstartilde)
	}
	if(reference.method == "population" & !maternal) {
		.ut <- cbind(.p, .Rtilde, .Rstartilde)	#
	}
	if(reference.method == "reciprocal" & !maternal) {
		.ut <- cbind(.p, .Fi, .Fi.startilde)	#
	}
	if(reference.method == "ref.cat" & maternal) {
		.ut <- cbind(.p, .R, .Rstartilde, .Rm, .Rmstartilde)
	}
	if(reference.method == "population" & maternal) {
		.ut <- cbind(.p, .Rtilde, .Rstartilde, .Rmtilde, .Rmstartilde)	#
	}	
	if(reference.method == "reciprocal" & maternal) {
		.ut <- cbind(.p, .Fi, .Fi.startilde, .Fim, .Fim.startilde)	#
	}

	if(!maternal) dimnames(.ut)[[2]] <- c(paste("p", 1:n.all, sep = ""), paste("RRc", 1:n.all, sep = ""), paste("RRcdd", 1:n.all, sep = ""))
	else
	dimnames(.ut)[[2]] <- c(paste("p", 1:n.all, sep = ""), paste("RRc", 1:n.all, sep = ""), paste("RRcdd", 1:n.all, sep = ""), paste("RRm", 1:n.all, sep = ""), paste("RRmdd", 1:n.all, sep = ""))

	return(.ut)
}
