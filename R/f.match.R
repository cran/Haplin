f.match <- function(d1, d2){

	if(!identical(names(d1), names(d2))) stop()
if(F){
		f.vis(system.time({
	.row.no <- seq(length.out = nrow(d2))
	.t <- tapply(.row.no, d2, identity)

	.mat1 <- as.matrix(d1)
	mode(.mat1) <- "character"

	.pos <- .t[.mat1]
	.pos1 <- sapply(.pos, function(x)x[1])
}))
}
.tag1 <- f.create.tag(d1, sep = "-")
.tag2 <- f.create.tag(d2, sep = "-")
.pos <- match(.tag1, .tag2)

if(F) if(!identical(.pos, .pos1)) stop()

return(.pos)

}