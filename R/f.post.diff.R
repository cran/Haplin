f.post.diff <- function(coeff, covar){
##
## FOR THE HAPLOTYPE FREQUENCIES, SUBTRACT FIRST PARAMETER FROM THE REST,
## TO "NORMALIZE". DO THIS BY MULTIPLYING WITH DIFFERENCE MATRIX
.names <- names(coeff[[1]])
.mf <- grep("mf", .names, value = F)

.D <- diag(length(coeff[[1]]))
dimnames(.D) <- list(.names, .names)
#.D[1:length(.mf), 1] <- -1
.D[.mf, .mf[1]] <- -1
#.D <- .D[-1,]
.D <- .D[-.mf[1],]
#
.vis <- F
f.vis(.D, vis = .vis)
f.vis(coeff[[1]], vis = .vis)
f.vis(.D %*% coeff[[1]], vis = .vis)
#
f.vis(.coef <- lapply(coeff, function(x) .D %*% x), vis = .vis)
f.vis(covar, vis = .vis)
f.vis(.cov <- lapply(covar, function(x) .D %*% x %*% t(.D)), vis = .vis)



return(list(coeff = .coef, covar = .cov))


}
