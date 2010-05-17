f.read.data <- function(indata, sep = " ", allele.sep=";", na.strings="NA", use.missing = F, markers="ALL", variables = 0, family = "mfc"){
## INDATA  : NAME OF ASCII FILE WITH ALLELE DATA
## SEP : SEPARATOR BETWEEN COLUMNS IN INDATA
## ALLELE.SEP : SEPARATOR WITHIN COMLUMNS
## NA.STRINGS : MISSING VALUES
## USE.MISSING : ENABLES YOU TO INCLUDE MISSING VALUES
## MARKERS : ENABLES YOU TO ONLY USE SOME LOCI OR ONLY ONE LOCUS
## VARIABLES : GIVES NUMBER OF VARIABLES (COVARIATES, CASE/CONTROL ETC) COLUMNS
## FAMILY : CHARACTER VARIABLE GIVING FAMILYDESIGN. "MFC" GIVES TRIAD DATA, "C" GIVES CASE/CONTROL DATA
## NB! NA.STRINGS CAN'T BE CODED EQUAL TO VALUES IN DATASET
## IF THE SEPARATOR IS SOME SORT OF WHITE SPACE, THE SEP ARGUMENT IN COUNT.FIELDS AND SCAN ARE LEFT OUT:
#
## PREPARE
if(allele.sep==".") allele.sep <- "\\."
if(na.strings == sep & sep == allele.sep) stop("You can't use equal coding for na.strings, sep and allele.sep!!")
#
##
if(na.strings == " " || na.strings == ""){
	## COUNT FIELDS, SET UP NAMES:
	.count <- count.fields(file = indata,sep = sep)
	.nvar <- .count[1]
	if(any(.count != .nvar)) stop("Number of variables differ from row to row")
	.data <- read.table(file = indata,sep=sep,na.strings=na.strings)
} else{
	if(sep == allele.sep){
		.count <- count.fields(file = indata,sep = sep)
		.nvar <- .count[1]
		if(any(.count != .nvar)){
			.mid.data <- scan(file = indata,sep = sep)
			stop("Sorry, something's wrong, missing function index.rowcol")
			.index <- NULL # BARE FOR AA UNNGAA FEILMELDING FRA R CMD check
			###.index <- index.rowcol(.mid.data,.mid.data == na.strings,which = "rows")
			if(length(.index)==0) stop("Number of variables differ from row to row")
			.stop <- cumsum(.count)
			.start <- c(1,.stop[-(length(.stop))]+1)
			.data <- matrix(NA,ncol = max(.count),nrow = length(.count))
			.j <- 1
			for(i in 1:(length(.count))){
				if(.j <= length(.index) & .index[.j] >= .start[i] & .index[.j] <= .stop[i]){
					.data[i,] <- c(.mid.data[.start[i]:.index[.j]],na.strings,.mid.data[(.index[.j]+1):.stop[i]])
					.j <- .j + 1
				}else{
					.data[i,] <- .mid.data[.start[i]:.stop[i]]
				}
			}
        }else{
			.data <- read.table(file=indata,sep=sep,na.strings=na.strings)
        }
		.data <- apply(.data,2,as.character)
		.data <- replace(.data,.data==na.strings,NA)
		.length <- 0
		.xdata <- NULL
		if(variables >0) .xdata <- .data[,1:variables]
		.gen.data <- .data[,(1+variables):(dim(.data)[2])]
		.teller1 <- seq(1,((dim(.gen.data)[2])-1),by=2)
		.teller2 <- seq(2,(dim(.gen.data)[2]),by=2)
		.gen.data <- lapply(1:length(.teller1),function(i,x,allele.sep,.teller1,.teller2) paste(x[,.teller1[i]],allele.sep,x[,.teller2[i]],sep=""),x = .gen.data,allele.sep = allele.sep,.teller1 = .teller1,.teller2 = .teller2)
		if(variables > 0){
			.data <- data.frame(.xdata,.gen.data)
		}else{
		  .data <- data.frame(.gen.data)
		}
	}else{
		## COUNT FIELDS, SET UP NAMES:
		.count <- count.fields(file = indata)
		.nvar <- .count[1]
		if(any(.count != .nvar)) stop("Number of variables differ from row to row")
		if(sep == " " | sep == ""){
		  .data <- read.table(file = indata,na.strings=na.strings)
		}else{
		  .data <- read.table(file = indata,na.strings=na.strings,sep=sep)
		}
	}
}
##
##
##
.xdata <- NULL
.data <- apply(.data,2,as.character)
.data <- replace(.data, .data==paste(na.strings,allele.sep,na.strings,sep=""), NA)
.data <- replace(.data, is.na(.data), "NA")
.data.full <- .data
row.names(.data.full) <- 1:(dim(.data.full)[1])
row.names(.data) <- 1:dim(.data)[1]
na.strings <- "NA"
.data <- apply(.data,2,as.character)
if(variables > 0){
	## SELECT VARIABLES COLUMNS
	###.xdata <- apply(.data[,1:variables,drop=F],2,as.character)
	.xdata <- .data[, 1:variables, drop=F]
	.xnamevec <- paste("x", 1:variables, sep = "")
	###.xnamevec <- c(sapply(1:variables,function(i) c(paste("x",i,sep=""))))
}
## GENETIC DATA:
.data.gen <- .data[,(1 + variables):dim(.data)[2], drop = F]
## TEST IF DATA CONTAIN MISSING VALUES:
.true <- apply(.data.gen,2,function(x,na.s,na.sep) x==paste(na.s,na.s,sep=na.sep)|is.na(x)|x==na.s,na.s=na.strings,na.sep=allele.sep)
###.data.gen <- apply(.data.gen,2,as.character)
.data.gen <- replace(.data.gen,.data.gen==na.strings,paste(na.strings,na.strings,sep=allele.sep))
## IF THE LENGTH OF THE CHARACTER IS GREATER THAN 2:
if(all(nchar(.data.gen[.true == F]) > 2)){
	.data.gen <-lapply(1:dim(.data.gen)[2],function(i,data,allele.sep) do.call("rbind",strsplit(unlist(data[,i]),split=allele.sep)),data=.data.gen,allele.sep=allele.sep)
	.data.gen <- do.call("cbind",.data.gen)
}
else if(all(nchar(.data.gen[.true == F]) == 2)){
	.data.gen <- replace(.data.gen,.data.gen == paste(na.strings,allele.sep,na.strings,sep=""),paste(na.strings,na.strings,sep=""))
	## SPLIT ALL CHARACTER STRINGS IN HALF, PUT SIDE BY SIDE IN SEPARATE COLUMNS
	.dim.tmp <- dim(.data.gen)
	.nchar <- nchar(.data.gen)
	.data.gen <- rbind(substring(.data.gen, first = 1, last = .nchar/2), substring(.data.gen, first = .nchar/2 + 1, last = .nchar))
	dim(.data.gen) <- c(.dim.tmp[1], 2 * .dim.tmp[2])
}
## FIND THE NAMES OF THE UNIQUE ALLELES IN EACH LOCUS:
.data.gen <- apply(.data.gen,2,as.character)
## IF ARGUMENT FAMILY IS EQUAL TO "MFC" THEN THERE IS 6 COLUMNS FOR EACH LOCI, ELSE IF FAMILY IS
## EQUAL TO "C" THEN THER IS 2 COLUMNS FOR EACH LOCI
if(family == "mfc"){
	.t <- 6
}else{
	if(family != "c") stop(" Problem with family argument!")
	.t <- 2
}
## FIND THE NUMBER OF LOCI
.nloci <- dim(.data.gen)[2]/.t
#
## CHECK THAT markers ARE WITHIN RANGE OF DATA
if(is.numeric(markers)){
	if(max(markers) > .nloci) stop(paste('Only ', .nloci, ' locus/loci found in file, "markers" argument specifies a marker number not in file.', sep = ""))
}
if(is.numeric(markers)){
	## SELECT THE MARKERS SPECIFIED IN markers ARGUMENT 
	.sel <- as.numeric(t(outer((markers-1)*.t, 0:(.t - 1), "+")) + 1)
	.data.gen <- .data.gen[, .sel, drop = F]
	###.tmp <- lapply(1:length(markers), function(i,markers,.data.gen,.t) .data.gen[,(1+(markers[i]-1)*.t):(markers[i]*.t)],markers=markers,.data.gen=.data.gen,.t=.t)
	###.data.gen <- do.call("cbind", .tmp)
	.nloci <- length(markers)
}
row.names(.data.gen) <- 1:(dim(.data.gen))[1]
#
## REMOVE ROWS WITH MISSING:
if(!use.missing){
	## CAN'T USE NA.EXCLUDE DUE TO DATA ON CHARACTER FORMAT
	if(variables > 0) .xdata <- .xdata[apply(.data.gen,1,function(x,na.strings) all(x != na.strings),na.strings=na.strings),,drop = F]
	.data.gen <- .data.gen[apply(.data.gen,1,function(x,na.strings) all(x != na.strings),na.strings=na.strings),]
	.na.message <- dim(.data.full)[1]-dim(.data.gen)[1]
}else{
	.tmp <- ((.data.gen == "NA") + 0) %*% rep(1,dim(.data.gen)[2]) == 0 # Q & D
	.na.message <- dim(.data.full)[1]-sum(.tmp)
}
.row.names <- row.names(.data.gen)
## FIND HOW MANY ROWS REMOVED DUE TO MISSING VALUES:
.omitted <- seq(nrow(.data.full))[is.na(match(row.names(.data.full),row.names(.data.gen)))]
if(length(.omitted) == 0) .omitted <- NULL
#
## CONSTRUCT A NAME VECTOR:
.namevector <- NULL
if(is.numeric(markers)){
	if(family == "mfc"){
		.namevector <- c(sapply(1:.nloci,function(i,markers) c(paste("l",markers[i],".m1",sep=""),paste("l",markers[i],".m2",sep=""),paste("l",markers[i],".f1",sep=""),paste("l",markers[i],".f2",sep=""),paste("l",markers[i],".c1",sep=""),paste("l",markers[i],".c2",sep="")),markers=markers))
	}else{
		.namevector <- c(sapply(1:.nloci,function(i,markers) c(paste("l",markers[i],".c1",sep=""),paste("l",markers[i],".c2",sep="")),markers=markers))
	}
}else{
	if(family == "mfc"){
	.namevector <- c(sapply(1:.nloci,function(i) c(paste("l",i,".m1",sep=""),paste("l",i,".m2",sep=""),paste("l",i,".f1",sep=""),paste("l",i,".f2",sep=""),paste("l",i,".c1",sep=""),paste("l",i,".c2",sep=""))))
	}else{
		.namevector <- c(sapply(1:.nloci,function(i) c(paste("l",i,".c1",sep=""),paste("l",i,".c2",sep=""))))
	}
}
## IF THERE ARE ANY VARIABLES IN DATA, THEY ARE ADDED :
if(variables > 0){
	.namevector <- c(.xnamevec, .namevector)
	.data.out <- cbind(apply(.xdata[,,drop=F],2,as.character), .data.gen)
}else .data.out <- .data.gen
if(dim(.data.out)[2] != length(.namevector)){
	## SOMEWHAT AD-HOC TEST FOR A MISSPECIFIED n.vars ARGUMENT
	cat("\n")
	stop("There's a problem with the number of variables in the data file.\n Are you sure the argument 'n.vars' is set to the correct value?", call. = F)
}
dimnames(.data.out) <- list(NULL, .namevector)
#
## RETURNS NUMBER OF ROWS DROPPED:
attr(.data.out,"rows.with.na")  <- .na.message
## RETURNS WHICH ROWS DROPPED:
attr(.data.out,"rows.dropped")  <- as.numeric(.omitted)
if(identical(markers, "ALL")){
	attr(.data.out,"markers") <- 1:.nloci
}else attr(.data.out, "markers") <- markers
return(.data.out)
}
 
