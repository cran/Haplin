print.haplinStrat <- function(x, ...){
##
##
##
cat('List of length ', length(x), ' containing results from a stratified haplin analysis.\nUse "summary" function to get detailed information.\n', sep = '') 

for(i in seq(along = x)){
	cat("\n#### Stratum = ", names(x)[i], " ####\n", sep = "")
	print(x[[i]])
}
return(invisible(x))

}
