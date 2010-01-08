print.info <- function(info){

.f.printlist <- function(x){
	.ut <- rep(NA, length(x))
	for(i in seq(along = x)){
		.ut[i] <- paste(' ', names(x)[i], ': "',  paste(x[i], collapse = '" "'), '"', sep = '')
	}
	return(.ut)
}


cat("\nHaplin call summary:\n")
cat('Filename: "', info$filename, '"\n', sep = "")
cat("File specifications:\n")
cat(.f.printlist(info$filespecs), sep = "\n")
cat("Model specifications:\n")
cat(.f.printlist(info$model), sep = "\n")
cat("Variables specifications:\n")
cat(.f.printlist(info$variables), sep = "\n")
cat("\n")




}
