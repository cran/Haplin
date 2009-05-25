f.prep.data.new <-  function(data, short = F)
  {
## MAKE SURE IT'S CHARACTER DATA, REGARDLESS:
.data <- apply(data,2,as.character)

## TEST WHETHER FAMILY IS "MFC" OR "C" (SHOULD RATHER HAVE DESIGN ARGUMENT)
.family <- any(regexpr("\\.m",dimnames(.data)[[2]])>0)

if(is.null(attr(data,"markers")))
      {
        markers <- "ALL"
      }
else
      {
        markers <- attr(.data,"markers")
      }

if(.family)
      {
        .t <- 6
      }
    else
      {
        .t <- 2
      }

#
#
#
    variables <- length(dimnames(.data)[[2]][regexpr("\\x",dimnames(.data)[[2]])>0])

    if(variables > 0 )
      {
## CONVERT VARIABLE PART TO LIST:
	.xdata <- f.matrix.to.list(.data[,1:variables, drop = F])
## STRIP WHITE SPACE IN BEGINNING AND END, TO DETECT NAs CORRECTLY, ETC.:
	.xdata <- lapply(.xdata, f.strip.white)
        .xnamevec <- dimnames(.data)[[2]][regexpr("\\x",dimnames(.data)[[2]])>0]        
#
## RECODE VARIABLES DATA :
.freq.x <- lapply(.xdata,table, exclude = c("NA", NA))
names(.freq.x) <- .xnamevec
.unique.x <- lapply(.freq.x,names)
.length.x <- as.vector(unlist(lapply(.unique.x,length)))
###        if(variables > 1) cat("Problem: Too many covariates!","\n")
###        if(any(.length.x > 2)) cat("Problem: Too many levels of covariate(s)!","\n")
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
    
    if(is.numeric(markers))
      {
        .alleles <- .frequency[markers]
      }
    else
      {
        .alleles <- .frequency
      }


	if(short) {

		if(variables > 0){
			.data <- cbind(.xdata, .data)
			attr(.data,"variables") <- .freq.x
			dimnames(.data)[[2]][1:variables] <- .xnamevec
		}			

			
	
	attr(.data, "alleles") <- .alleles
	    attr(.data,"rows.with.na")  <- attr(data,"rows.with.na")
	    	## RETURNS WHICH ROWS DROPPED:
	    attr(.data,"rows.dropped")  <- attr(data,"rows.dropped")
	return(.data)
}










########
########
## TROR IKKE JEG TRENGER NOE AV DET FOLGENDE.




	.nrows <- dim(.data)[1]

	.mendel.pos.list <- vector(.nrows, mode = "list")
	.mendelian.consistency <- vector(.nrows, mode = "logical")
	
f.vis(.nallele, vis = T)

	for (i in 1:.nrows){
		if(i%%10 == 0) cat(i, " ")
		.mendel.pos <- vector(length = .nloci, mode = "list")
		for (j in 1:.nloci){
			.mid <-  .data[i,(1+(j-1)*.t):(.t*j), drop = F]
			.mendel.pos[[j]] <- f.sort.alleles.new(.mid, A = length(.alleles[[j]]))
		}
		.mendelian.consistency[i] <- all(sapply(.mendel.pos, function(x) attr(x, "mendelian.consistency")))
		.mendel.pos.list[[i]] <- .mendel.pos
	}

	cat("\n")

##	.mendelian.consistency <- lapply(.mendelian.consistency, function(x)
	
	
##	apply(do.call("cbind", .mendelian.consistency), 1, all)

##	return(.mendelian.consistency)

###	.mendel.pos.list <- .mendel.pos.list[.mendelian.consistency]


## return(.mendel.pos.list)


##	.vekt <- cumprod(.nallele ^4)
##	.vekt <- c(1, .vekt[-length(.vekt)])
	
##	return(.mendel.pos.list)
	
if(F)	.f.compute.pos <- function(pos.list, vekt){
		for (i in seq(along = pos.list)){
			pos.list[[i]] <- (pos.list[[i]] - 1) * vekt[i]
		}
		pos.list <- do.call("expand.grid", pos.list)
		pos.list <- apply(pos.list, 1, sum)
		pos.list + 1
	}	
		
	.f.compute.pos <- function(pos.list, nallele){
		.pos.list.tmp <- vector(length(pos.list), mode = "list")
		pos.list <- do.call("expand.grid", pos.list)
		for (i in seq(along = pos.list)){
			f.vis(pos.list[,i])
			f.vis(nallele[i])
			.pos.list.tmp[[i]] <- f.pos.in.grid(pos = pos.list[,i], A = rep(nallele[i], 4))
		}

		.ut <- do.call("cbind", .pos.list.tmp)		



##		pos.list <- apply(pos.list, 1, sum)
##		pos.list + 1
		return(.ut)
	}	


##	.mendel.pos.list <- lapply(.mendel.pos.list, .f.compute.pos, vekt = .vekt)
	.mendel.pos.list <- lapply(.mendel.pos.list, .f.compute.pos, nallele = .nallele)


## return(.mendel.pos.list)


### return(.mendel.pos.list)
	.ind <- rep(seq(along = .mendel.pos.list), sapply(.mendel.pos.list, function(x) dim(x)[1]))
## return(.ind)	

	.haplocomb <- do.call("rbind", .mendel.pos.list)

## return(cbind(.ind, .haplocomb))
	
	.m1 <- .m2 <- .f1 <- .f2 <- numeric(length = length(.ind))
	
	.vekt <- cumprod(.nallele)
	.vekt <- c(1, .vekt[-length(.vekt)])
	
	f.vis(.vekt, vis = T)
	
	f.vis(.haplocomb[1:20,], vis = T)
	
	for (i in seq(along = .nallele)){
		f.vis(4*(i-1) + 1, vis = T)
		.m1 <- (.haplocomb[, 4*(i-1) + 1] - 1) * .vekt[i] + .m1
		.m2 <- (.haplocomb[, 4*(i-1) + 2] - 1) * .vekt[i] + .m2
		.f1 <- (.haplocomb[, 4*(i-1) + 3] - 1) * .vekt[i] + .f1
		.f2 <- (.haplocomb[, 4*(i-1) + 4] - 1) * .vekt[i] + .f2
	}

##	.haplocomb <- f.pos.to.haplocomb(unlist(.mendel.pos.list), prod(.nallele)) # GAAR OVER BEKKEN ETTER VANN, MEN....
	
	.ut <- as.data.frame(cbind(m1 = .m1 + 1, m2 = .m2 + 1, f1 = .f1 + 1, f2 = .f2 + 1, ind = .ind))

	attr(.ut, "alleles") <- .alleles
	attr(.ut,"rows.with.Mendelian.inconsistency")  <- which(!.mendelian.consistency)

return(.ut)

        for(j in 1:.nloci)
          {
	## PICKS OUT DATA FOR LOCUS J:
            .mid <-  .data[,(1+(j-1)*.t):(.t*j), drop = F]
		.mendel.sort <- vector(length = dim(.mid)[1], mode = "list")
		for (i in seq(length = dim(.mid)[1]))
		{
			.mendel.sort[[i]] <- f.sort.alleles.new(.mid[i,], A = length(.alleles[[j]]))
		}
			.mendelian.consistency[[j]] <- sapply(.mendel.sort, function(x) attr(x, "mendelian.consistency"))
			.mendel.sort.list[[j]] <- .mendel.sort
		}		










	.ut <- lapply(.mendel.sort.list, function(x, i) x[i], .mendelian.consistency)



	










	attr(.ut, "alleles") <- .alleles
	attr(.ut,"rows.with.Mendelian.inconsistency")  <- which(!.mendelian.consistency)









	return(.ut)


    return(.data)




#################################
#################################





    
    ## SORT ALLELES ROW BY ROW AND PUT INTO A DATA FRAME:
    .nydata <- NULL
    .any <- NULL
    .row.names <- 1:(dim(.data)[1])
    if(.family)
      {
        .namevector <- dimnames(.data)[[2]][regexpr("\\.c",dimnames(.data)[[2]])<0]
        for(j in 1:.nloci)
          {
            ## PICS OUT DATA FOR LOCUS J:
            .mid <-  .data[,(1+(j-1)*.t):(.t*j)]
            .mid <- t(apply(.mid,1,f.sort.alleles))
            dimnames(.mid) <- list(.row.names,NULL)
            ## ROW NUMBER RETURNED IF F.SORT.ALLELES FIND ANY ROWS WITH MENDELIAN INCONSISTENCY
            .mendelian <- do.call("c",lapply(1:(dim(.mid)[1]),function(i,x) attr(f.sort.alleles(x[i,]),"mendelian.consistency"),x=.mid))
            .mid.mend <- .mid[.mendelian,1:4]
            .mid <- .mid[,1:4]
            dimnames(.mid) <- list(.row.names,NULL)
            ## FIND IF THERE ARE ROWS WITH MENDELIAN INCONSISTENCIES
            .any.mend <- match(dimnames(.mid)[[1]],dimnames(.mid.mend)[[1]])
            .any <- c(.any,dimnames(.mid)[[1]][is.na(.any.mend)])
            .nydata <- cbind(.nydata,as.matrix(.mid))
          }
      }
    else
      {
        .namevector <- dimnames(.data)[[2]]
        .nydata <- .data
        .mid <- .data
      }
    dimnames(.nydata)[[2]] <- .namevector
    .match <- match(dimnames(.mid)[[1]],.any)
    .any.rm <- (1:length(.match))[!is.na(.match)]
    if(length(.any.rm)>0)
      {
        .nydata <- .nydata[-(.any.rm),]
      }
    ## COMPUTE THE FREQUENCIES OF THE DIFFERENT ALLELE CONFIGURATIONS:
    ## MAKE SURE THAT .NYDATA IS NOT A FACTOR
    .nydata <- apply(.nydata,2,as.numeric)
    .nydata <- data.frame(.nydata)
    if(variables > 0)
      {
        if(length(.any.rm) > 0) .xdata <- .xdata[-(.any.rm),]
        .namevector <- c(.xnamevec,.namevector)
        .nydata <- cbind(.xdata,.nydata)
      }
    .nydata <- apply(.nydata,2,as.numeric)
    dimnames(.nydata) <- list(NULL,.namevector)
    
##    return(markers)


    
    
    return(.nydata)
    
    
    
    
    ## FINDS THE FREQUENCIES OF THE ALLELE COMBINATIONS:
    .aggregate <- apply(.nydata,1,function(x) paste(x,collapse="&"))
    .freq <- as.vector(table(.aggregate))
    .aggregate <- names(table(.aggregate))
    .aggregate <- do.call("cbind",unpaste(.aggregate,sep="&"))
    dimnames(.aggregate) <- list(NULL,.namevector)
    .aggregate <- cbind(.aggregate[,,drop=F],freq = .freq)
    .aggregate <- apply(.aggregate,2,as.numeric)
    if(is.null(dim(.aggregate)))
      {
        .data <- t(as.matrix(.aggregate))
      }
    else
      {
        .data <- .aggregate
      }
    if(!is.null(dim(.data)))
      {
        .na <- apply(.data,1,function(x) any(is.na(x)))
        .data <- .data[order(.na),,drop=F]
      }
    ## RETURNS ALLELENAMES AS AN ATTRIBUTE TO .DATA:
    if(is.numeric(markers))
      {
        attr(.data,"alleles") <- .frequency[markers]
      }
    else
      {
        attr(.data,"alleles") <- .frequency
      }
    ## RETURNS LINES WITH MED MENDELIAN INCON.
    attr(.data,"rows.with.Mendelian.inconsistency")  <- as.numeric(.any)
    ## RETRUNS DETAILS OF THE COVARIATES, IF COVARIATES:
    if(variables > 0)
      {
        names(.freq.x) <- .xnamevec
        attr(.data,"variables") <- .freq.x
      }
    attr(.data,"rows.with.na")  <- attr(data,"rows.with.na")
    ## RETURNS WHICH ROWS DROPPED:
    attr(.data,"rows.dropped")  <- attr(data,"rows.dropped")
    if(!is.null(attr(data,"markers"))) attr(.data,"markers") <- attr(data,"markers")
    return(.data)
  }
