f.strip.white <- function(x, R = T){
##
## REMOVES WHITE SPACE (SPACE AND TABS) AT BEGINNING AND END OF EACH STRING
## ELEMENT OF THE VECTOR x
##

## KANSKJE LIKE GREIT AA BRUKE gsub (I R), MEN DA FORSVINNER
## OGSAA INNMAT SOM ER WHITE.
##

	if(R) .x.split <- strsplit(x, split = "")
	else {
		.x.split <- lapply(x, function(x){
				substring(x, first = 1:nchar(x), last = 1:nchar(x))
			})
	}
		
	
	.f.find.white <- function(x){
	## x IS A VECTOR, RESULT IS A LOGICAL, TRUE WHERE THERE IS WHITE SPACE	
				
		.white <- (x == " ") | (x == "\t")
		
		if(all(.white)) return("") # ALSO COVERS EMPTY STRINGS

		.which.F <- which(!.white)
		.range <- range(.which.F)
		
		return(paste(x[.range[1]:.range[2]], collapse = ""))		
		
	}
	
	
	.res <- sapply(.x.split, .f.find.white)
	
	return(.res)


}
