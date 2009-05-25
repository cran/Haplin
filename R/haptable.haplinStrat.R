haptable.haplinStrat <- function(object){

.tab <- lapply(object, haptable)

for(i in seq(along = .tab)){
	.tab[[i]] <- cbind(stratum = names(.tab)[i], .tab[[i]])
}


return(.tab)

}
