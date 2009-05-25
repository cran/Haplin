f.post.diff <- function(coeff, covar){
#
## FOR THE HAPLOTYPE FREQUENCIES, SUBTRACT FIRST PARAMETER FROM THE REST,
## TO "NORMALIZE". DO THIS BY MULTIPLYING WITH DIFFERENCE MATRIX
.names <- names(coeff[[1]])
.mf <- grep("mf", .names, value = T)

.D <- diag(length(coeff[[1]]))
dimnames(.D) <- list(.names, .names)
.D[1:length(.mf), 1] <- -1
.D <- .D[-1,]
#
f.vis(.D)
f.vis(coeff[[1]])
f.vis(.D %*% coeff[[1]])
#
f.vis(.coef <- lapply(coeff, function(x) .D %*% x))
f.vis(covar)
f.vis(.cov <- lapply(covar, function(x) .D %*% x %*% t(.D)))


return(list(coeff = .coef, covar = .cov))


}