"f.posttest"<-
function(coef, cov, mf, c, cdd, m, mdd, test)
{
.vis <- F
.coef <- coef
.cov <- cov
## SELECT PARAMETERS TO BE TESTED
## test <- c("child.single", "child.double", "maternal.single", "maternal.double")

.names.list <- list(haplo.freq = mf, child = c(c, cdd), maternal = c(m, mdd))
.nam <- names(.names.list)
if(!all(test %in% .nam)) stop('Invalid input in argument "test"')
f.vis(.velg <- unlist(.names.list[.nam %in% test]), vis = .vis) # MAKE SURE SELECTION IS IN CORRECT ORDER


#
f.vis(.coef <- lapply(.coef, function(x)x[.velg, , drop = F]), vis = .vis)
f.vis(.cov <- lapply(.cov, function(x) x[.velg, .velg, drop = F]), vis = .vis)
#
## RESHAPE COEFFICIENTS AND COVARIANCE MATR. INTO FULL SIZE
.n.pars <- length(.coef[[1]])
.l <- length(.coef)
f.vis(.coef.vec <- unlist(.coef), vis = .vis)
f.vis(.cov.mat <- f.bdiag(.cov), vis = .vis)
#
## BUILD CONTRAST MATRIX
.A <- f.post.contrasts(test.type = "interaction", n.res = .l, n.pars = .n.pars)

#
## DO CHI-SQUARED TEST
.chisq.res <- f.post.chisq(coeff = .coef.vec, covar = .cov.mat, contrast.mat = .A)


return(.chisq.res)

}
