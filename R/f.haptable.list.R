f.haptable.list <- function(object){
##
## APPLIES haptable TO EACH ELEMENT OF object
## object IS A list OF haplin OBJECTS
## 

.tab <- lapply(object, haptable)

for(i in seq(along = .tab)){
	.tab[[i]] <- cbind(stratum = names(.tab)[i], .tab[[i]])
}

.tab <- do.call("rbind", .tab)
rownames(.tab) <- NULL

return(.tab)

}
