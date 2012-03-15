f.prep.data <-  function(data, info)
  {
## MAKE SURE IT'S CHARACTER DATA, REGARDLESS:
.data <- apply(data,2,as.character)
#
## TEST WHETHER FAMILY IS "MFC" OR "C" (USING 6 VERSUS 2 COLUMNS)
if((info$model$design == "triad") | (info$model$design == "cc.triad")) {
	.t <- 6
}else if(info$model$design == "cc") .t <- 2
#
##
markers <- info$filespecs$markers
variables <- info$filespecs$n.vars
#
if(variables > 0){
	## CONVERT VARIABLE PART TO LIST:
	.xdata <- f.matrix.to.list(.data[,1:variables, drop = F])
	## STRIP WHITE SPACE IN BEGINNING AND END, TO DETECT NAs CORRECTLY, ETC.:
	.xdata <- lapply(.xdata, f.strip.white)
	.xnamevec <- dimnames(.data)[[2]][seq(length.out = variables)]
	#
	## RECODE VARIABLES DATA :
	.freq.x <- lapply(.xdata,table, exclude = c("NA", NA))
	names(.freq.x) <- .xnamevec
	.unique.x <- lapply(.freq.x,names)
	.length.x <- as.vector(unlist(lapply(.unique.x,length)))
	.xdata.ny <- .xdata
	.xdata.ny <- do.call("cbind",lapply(1:length(.unique.x),function(i,.xdata.ny,.xdata,.unique.x) f.recode(.xdata.ny[[i]],.xdata[[i]],.unique.x[[i]]),.xdata.ny=.xdata.ny,.xdata=.xdata,.unique.x=.unique.x))
	.xdata.ny[.xdata.ny == "NA"] <- NA
	.xdata <- apply(.xdata.ny,2,as.numeric)
}
#
#
## GENETIC DATA :
.data <- .data[,(1 + variables):(dim(.data)[2])]
.nloci <- (dim(.data)[2])/.t
.unique <- lapply(1:.nloci,function(i,data,.t) names(table(unlist(data[,(1+(i-1)*.t):(.t*i)]),exclude=c("NA",NA))),data=.data,.t=.t)
## FIND THE NUMBER OF ALLELES IN EACH LOCUS, THE LENGTH OF .NALLELE IS EQUAL TO THE NUMBER OF LOCI:
.nallele <- as.vector(unlist(lapply(1:.nloci,function(i,x) length(x[[i]]),x=.unique)))
## FIND THE FREQUENCY OF EACH ALLELE:
.sum <- lapply(1:.nloci,function(i,data,.t) table(data[,(1+(i-1)*.t):(.t*i)],exclude=c("NA",NA)),data=.data,.t)

## LIST OF FREQUENCY OF EACH ALLELE, LIST HAS LENGTH EQUAL TO THE NUMBER OF LOCI:
.frequency <- lapply(1:length(.sum),function(i,sum,unik) f.names(i,sum,unik) ,sum=.sum,unik=.unique)
## RECODE GENETIC DATA :
.data.ny <- .data
.data.ny <- do.call("cbind",lapply(1:.nloci,function(i,.data.ny,.data,.unique,.t) f.recode(.data.ny[,(1+(i-1)*.t):(.t*i)],.data[,(1+(i-1)*.t):(.t*i)],.unique[[i]]),.data.ny=.data.ny,.data=.data,.unique=.unique,.t=.t))
.data.ny[.data.ny == "NA"] <- NA
.data <- apply(.data.ny,2,as.numeric)
#
if(variables > 0){
	.data <- cbind(.xdata, .data)
	attr(.data,"variables") <- .freq.x
	dimnames(.data)[[2]][1:variables] <- .xnamevec
}			
#
names(.frequency) <- markers
attr(.data, "alleles") <- .frequency
attr(.data,"rows.with.na")  <- attr(data,"rows.with.na")
## RETURNS WHICH ROWS DROPPED:
attr(.data,"rows.dropped")  <- attr(data,"rows.dropped")
return(.data)

}
