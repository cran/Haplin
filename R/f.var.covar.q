f.var.covar <- function(pred, X, data, info){
##
## COMPUTES EXPLICITLY THE VARIANCE-COVARIANCE MATRIX FOR THE LIKELIHOOD WITH MISSING DATA
##
## X IS THE DESIGN MATRIX (COMPLETE GRID), pred ARE PREDICTED FREQUENCIES FROM LAST RUN OF THE EM,
## data ARE THE "ORIGINAL" DATA USED BY EM
##
## NOTE: SCORE IS COMPUTED ACCORDING TO MULTINOMIAL SAMPLING MODEL, THUS ONE OF THE FIRST COLUMNS IS REALLY REDUNDANT
## THE COMPUTATIONS OF VAR-COVAR FROM MULTINOMIAL MODEL NOT ALWAYS SUCCESSFUL, SHOULD AT LEAST BE DONE ON X[,-1]
#
#
##
.n.sel.haplos <- sum(info$haplos$selected.haplotypes)
.design <- info$model$design
.xchrom <- info$model$xchrom
#
## STANDARD POISSON CONTRIBUTION:
.pX <- pred * X
.d2Poisson <- -t(X) %*% (.pX) # THIS IS STANDARD POISSON WITH VALUES PREDICTED ACCORDING TO EM
##
##
#
## MIDLERTIDIG SJEKK: #-#
.sjekk <- tapply(data$orig.lines, data$ind, unique)
.sjekk0 <- table(.sjekk)
.sjekk1 <- table(tapply(data$ind, data$orig.lines, unique))
if(any(.sjekk0 != 1) | any(.sjekk1 != 1)) stop() # SJEKKER AT ind OG orig.lines ER EN-TIL-EN, BURDE IKKE VÆRE NØDVENDIG

#
## "SLOW" COMPUTATION, TRIAD BY TRIAD:
if(F){#
	.ind <- data$ind
	.unique.ind <- unique(.ind)
	.m <- dim(X)[1]		
	.matsum <- matrix(0, ncol = .m, nrow = .m)	
	for(i in seq(along = .unique.ind)){	
	cat(i,"\n")
		.sel <- .ind == .unique.ind[i]
###		.pos <- f.pos.in.grid(A = rep(.n.sel.haplos, 4), comb = data[.sel,1:4])
	stop("Tviler paa at de folgende linjene er rett!")
	if((.design == "triad") & !.xchrom){
		.pos <- f.pos.in.grid(A = rep(.n.sel.haplos, 4), comb = as.matrix(data[,c("m1", "m2", "f1", "f2")]))
	}
	if((.design == "triad") & .xchrom){
		.pos <- f.pos.in.grid(A = rep(.n.sel.haplos, 3), comb = as.matrix(data[,c("m1", "m2", "f2")]))
	}
	if(.design == "cc"){
		.pos <- f.pos.in.grid(A = c(rep(.n.sel.haplos, 2), 2), comb = as.matrix(data[.sel,c("c1", "c2", "cc")]))
	}
	if(.design == "cc.triad"){
		.pos <- f.pos.in.grid(A = c(rep(.n.sel.haplos, 4), 2), comb = as.matrix(data[.sel,c("m1", "m2", "f1", "f2", "cc")]))
	}

		.gamma <- rep(0, .m)
		.gamma[.pos] <- pred[.pos]
		.gamma <- .gamma/sum(.gamma)		
		.mat <- diag(.gamma) - .gamma %*% matrix(.gamma, nrow = 1)
		.matsum <- .matsum + .mat
	}
	.d2phi.del1.tmp <- t(X) %*% .matsum %*% X
}

##
##
#
## MATCH EVERYTHING TO data
if(F){# ERSTATTET AV f.pos.match
	if((.design == "triad") & !.xchrom){
		.pos <- f.pos.in.grid(A = rep(.n.sel.haplos, 4), comb = as.matrix(data[,c("m1", "m2", "f1", "f2")]))
	}
	if((.design == "triad") & .xchrom){
		.pos <- f.pos.in.grid(A = c(rep(.n.sel.haplos, 3), 2), comb = as.matrix(data[,c("m1", "m2", "f2", "sex")]))
	}
	if(.design == "cc"){
		if(.xchrom)stop("Not implemented")
		.pos <- f.pos.in.grid(A = c(rep(.n.sel.haplos, 2), 2), comb = as.matrix(data[,c("c1", "c2", "cc")]))
	}
	if(.design == "cc.triad"){
		if(.xchrom)stop("Not implemented")
		.pos <- f.pos.in.grid(A = c(rep(.n.sel.haplos, 4), 2), comb = as.matrix(data[,c("m1", "m2", "f1", "f2", "cc")]))
	}
	.pos.test <- f.pos.match(data = data, design = .design, xchrom = .xchrom, n.sel.haplos = .n.sel.haplos)
	if(!all.equal(.pos, .pos.test)) stop()
}
.pos <- f.pos.match(data = data, design = .design, xchrom = .xchrom, n.sel.haplos = .n.sel.haplos)
#
.X <- X[.pos,]
.m <- dim(.X)[1]
.k <- dim(.X)[2]
.l <- pred[.pos]
.colnames <- colnames(.X)
#	
## NORMALIZE PREDICTIONS OVER AMBIGUITY GROUPS:
.ind <- data$ind #-#
.orig.lines <- data$orig.lines
#.lstar.gugg <- .l/f.groupsum(.l, INDICES = .ind) # NORMALIZED #-#
.lstar <- .l/f.groupsum(.l, INDICES = .orig.lines) # NORMALIZED
#if(any(abs(.lstar.gugg - .lstar) > 1e-5)) stop() #-#

#.Xl.gugg <- .lstar.gugg * .X #-#
.Xl <- .lstar * .X
#	
## COMPUTE INDIVIDUAL SCORE PARTS:	
#.lTX.gugg <- tapply(as.numeric(.Xl), list(ind = rep(.ind, .k), col = rep(1:.k, each = .m)), sum) ## THIS IS (gamma_i)^T X, WITH ONE ROW FOR EACH i AND ONE COLUMN FOR EACH COLUMN OF X. EACH ROW IS THE INDIVIDUAL PART OF THE SCORE VECTORS	#-#
.lTX <- tapply(as.numeric(.Xl), list(orig.lines = rep(.orig.lines, .k), col = rep(1:.k, each = .m)), sum) ## THIS IS (gamma_i)^T X, WITH ONE ROW FOR EACH i AND ONE COLUMN FOR EACH COLUMN OF X. EACH ROW IS THE INDIVIDUAL PART OF THE SCORE VECTORS
colnames(.lTX) <- .colnames
###if(any(abs(.lTX0 - .lTX) > 1e-10)) stop()#-#

#.d2phi.1.gugg <- t(.lTX.gugg) %*% .lTX.gugg #-#
#.d2phi.2.gugg <- t(.X) %*% (.Xl.gugg) #-#
#.d2phi.del1.gugg <- - .d2phi.1.gugg + .d2phi.2.gugg #-#
#
## COMPUTE FIRST AND SECOND PART OF AMBIGUITY-PART OF SECOND DERIV. MULTINOMIAL LOGLIKE
.d2phi.1 <- t(.lTX) %*% .lTX
.d2phi.2 <- t(.X) %*% (.Xl)
.d2phi.del1 <- - .d2phi.1 + .d2phi.2
#	
## ADD STANDARD POISSON AND PART DUE TO AMBIGUITIES:
#.d2.loglike.Poisson.gugg <- .d2phi.del1.gugg + .d2Poisson #-#
.d2.loglike.Poisson <- .d2phi.del1 + .d2Poisson
#
## SCORE COMPUTATION (COMPUTED FROM MULTINOMIAL FORMULA t(gamma_i - gamma) %*% X
## COMMON ELEMENT FOR ALL SCORE VECTORS:
.gammaTX <- (t(.pX/sum(pred)) %*% rep(1, dim(.pX)[1]))[,1] ## LAST SUBSETTING REDUCES MATRIX TO NAMED VECTOR
.score <- t(t(.lTX) - .gammaTX) # SUBTRACT .gammaTX COLUMNWISE
#
##
if(F){
## COMPUTE VAR-COVAR ACCORDING TO MULTINOMIAL FORMULATION:
###	.s1 <- .d2Poisson/sum(pred)
	.s1 <- .d2Poisson ## THIS SHOULD BE n t(X) %*% diag(gamma) %*% X = t(X) %*% diag(pred) %*% X
	.s2 <- sum(pred) * .gammaTX %*% t(.gammaTX)
	.inf.multinom <- -(.d2phi.del1 + .s1 + .s2) ## INFORMATION MATRIX FROM MULTINOMIAL W/AMBIGUITIES
	.var.covar.multinom <- solve(.inf.multinom) ## DETTE GÅR IKKE ALLTID BRA!
}
#
## INVERT TO OBTAIN VAR-COVAR:	
#.var.covar.Poisson.gugg <- -solve(.d2.loglike.Poisson.gugg) #-#
.var.covar.Poisson <- -solve(.d2.loglike.Poisson)
#if(any(abs(.var.covar.Poisson - .var.covar.Poisson.gugg) > 1e-4)) stop() #-#
#
##
###	return(list(var.covar = .var.covar.Poisson, score = .score, inf = -.inf.multinom))
### return(list(var.covar = .var.covar.Poisson, score = .score, var.covar.multinom = .var.covar.multinom))
return(list(var.covar = .var.covar.Poisson, score = .score))
}
