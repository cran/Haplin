f.pos.match <- function(data, info){
##
## FOR data, COMPUTES THE POSITION OF EACH LINE IN A COMPLETE
## GRID AS USED IN THE GLM
#
## EXTRACT CHARACTERISTICS FROM THE DESIGN MATRIX
.char <- f.design.make(info = info, ret.characteristics = T)
#
if(info$model$design == "cc"){
	## MATCH TO CORRECT COLUMNS FOR CASE-CONTROL
	names(.char)[names(.char) == "m2"] <- "c1"
	names(.char)[names(.char) == "f2"] <- "c2"
}

if((info$model$design == "cc") & info$model$xchrom){
	#
	## Different coding in data and grid
	.char$cc <- .char$cc + 1
	#
	## Match data to grid
	if(info$model$comb.sex %in% c("males", "females")){
		.pos <- f.match(data[, c("c1", "c2", "cc")], .char)
	}
	if(info$model$comb.sex %in% c("single", "double")){
		.char$sex <- .char$sex + 1
		.pos <- f.match(data[, c("c1", "c2", "sex", "cc")], .char)
	}
}else{
	## MATCH PREDICTED PROBABILITIES (GRID) TO ORIGINAL DATA:
	.pos <- f.pos.in.grid(A = .char, comb = as.matrix(data[, names(.char)]))
}
#
##
return(.pos)
}
