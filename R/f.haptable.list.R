f.haptable.list <- function(object){
##
## APPLIES haptable TO EACH ELEMENT OF object
## object IS A list OF haplin OBJECTS
## 
#
## REMOVE MISSING (FAILED) RUNS
object[is.na(object)] <- NULL
#
## CREATE STANDARD haptable FOR EACH ELEMENT
.tab <- lapply(object, haptable)
#
## ADD STRATUM NAME TO TABLE
for(i in seq(along = .tab)){
	.tab[[i]] <- cbind(stratum = names(.tab)[i], .tab[[i]])
}
#
## JOIN ROWS FROM DIFFERENT STRATA
.tab <- do.call("rbind", .tab)
rownames(.tab) <- NULL
#
return(.tab)

}
