f.beta.haplo.freq.asymp <- function(haplo.freq){ 
	.haplo.freq <- haplo.freq/sum(haplo.freq)
	.sum <- 1/.haplo.freq[1]
	.beta <- log(.haplo.freq[-1]*.sum)
	return(c(0,.beta))
}