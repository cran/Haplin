f.post.contrasts <- function(test.type, n.res, n.pars){



if(test.type == "interaction"){


	f.vis(.id <- diag(n.pars))
	f.vis(.A1 <- matrix(rep(1, n.res - 1)))

	f.vis(.A2 <- -diag(n.res - 1))
	f.vis(cbind(.A1, .A2))
	f.vis(.A <- cbind(.A1, .A2) %x% .id)
	
}

return(.A)

}