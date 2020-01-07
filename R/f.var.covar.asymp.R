
## Computes variance-covariance
f.var.covar.asymp <- function(X, data, pred, ncells, norig, orig, info){
	# .p <- pred/sum(pred) # this is not needed since standardized within .f.b anyway
	.b <- f.b.asymp(data = data, p = pred, ncells = ncells, norig = norig, orig = orig, info = info)
	.mat <- f.mat.asymp(.b, w = pred, ncells=ncells)
	.d2l <- t(X) %*% (.mat - diag(pred)) %*% X
	.var.covar <- -solve(.d2l)
	return(.var.covar)
}
