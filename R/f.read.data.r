f.read.data <-
  function(indata,sep = " ",allele.sep=";",na.strings="NA",use.missing = F,markers="ALL",variables = 0,family = "mfc")
{
  ## INDATA  : ASCII FILE WITH ALLELE DATA
  ## SEP : SEOARATOR BETWEEN COLUMNS IN INDATA
  ## ALLELE.SEP : SEPARATOR WITHIN COMLUMNS
  ## NA.STRINGS : MISSING VALUES
  ## USE.MISSING : ENABLES YOU TO INCLUDE MISSING VALUES
  ## MARKERS : ENABLES YOU TO ONLY USE SOME LOCI OR ONLY ONE LOCUS
  ## VARIABLES : GIVES NUMBER OF VARIABLES (COVARIATES, CASE/CONTROL ETC) COLUMNS
  ## FAMILY : CHARACTER VARIABLE GIVING FAMILYDESIGN. "MFC" GIVES TRIADE DATA, "C" GIVES CASE/CONTROL DATA
  ## NB! NA.STRINGS CAN'T BE CODED EQUAL TO VALUES IN DATASET
  ## IF THE SEPARATOR IS SOME SORT OF WHITE SPACE, THE SEP ARGUMENT IN COUNT.FIELDS AND SCAN ARE LEFT OUT:


  if(allele.sep==".") allele.sep <- "\\."


  if(na.strings == sep & sep == allele.sep)
    {
      stop("You can't use equal coding for na.strings, sep and allele.sep!!")
    }


if(na.strings == " " || na.strings == "")
        {
          ## COUNT FIELDS, SET UP NAMES:
          .count <- count.fields(file = indata,sep = sep)
          .nvar <- .count[1]
          if(any(.count != .nvar))
            stop("Number of variables differ from row to row")
          .data <- read.table(file = indata,sep=sep,na.strings=na.strings)
        }
      else
        {
          if(sep == allele.sep)
            {
              .count <- count.fields(file = indata,sep = sep)
              .nvar <- .count[1]
              if(any(.count != .nvar))
                {
                  .mid.data <- scan(file = indata,sep = sep)
                  .index <- index.rowcol(.mid.data,.mid.data == na.strings,which = "rows")
                  if(length(.index)==0)
                    stop("Number of variables differ from row to row")
                  .stop <- cumsum(.count)
                  .start <- c(1,.stop[-(length(.stop))]+1)
                  .data <- matrix(NA,ncol = max(.count),nrow = length(.count))
                  .j <- 1
                  for(i in 1:(length(.count)))
                    {
                      if(.j <= length(.index) & .index[.j] >= .start[i] & .index[.j] <= .stop[i])
                        {
                          .data[i,] <- c(.mid.data[.start[i]:.index[.j]],na.strings,.mid.data[(.index[.j]+1):.stop[i]])
                          .j <- .j + 1
                        }
                      else
                        {
                          .data[i,] <- .mid.data[.start[i]:.stop[i]]
                        }
                    }
                }
              else
                {
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
              if(variables > 0)
                {
                  .data <- data.frame(.xdata,.gen.data)
                }
              else
                {
                  .data <- data.frame(.gen.data)
                }
            }
          else
            {
              ## COUNT FIELDS, SET UP NAMES:
              .count <- count.fields(file = indata)
              .nvar <- .count[1]
              if(any(.count != .nvar))
                stop("Number of variables differ from row to row")
              if(sep == " " | sep == "")
                {
                  .data <- read.table(file = indata,na.strings=na.strings)
                }
              else
                {
                  .data <- read.table(file = indata,na.strings=na.strings,sep=sep)
                }
            }
        }
  .xdata <- NULL
  .data <- apply(.data,2,as.character)
  .data <- replace(.data,.data==paste(na.strings,allele.sep,na.strings,sep=""),NA)
  .data <- replace(.data,is.na(.data),"NA")
  .full.data <- .data
  row.names(.full.data) <- 1:(dim(.full.data)[1])
  row.names(.data) <- 1:dim(.data)[1]
  na.strings <- "NA"
  .data <- apply(.data,2,as.character)
  if(variables > 0)
    {
      ## PICS OUT THE VARIABLES COLUMNS
      .xdata <- apply(.data[,1:variables,drop=F],2,as.character)
      .xnamevec <- c(sapply(1:variables,function(i) c(paste("x",i,sep=""))))
    }
  ## GENETIC DATA:
  .data <- .data[,(1 + variables):dim(.data)[2], drop = F]
  ## TEST IF DATA CONTAIN MISSING VALUES:
  .true <- apply(.data,2,function(x,na.s,na.sep) x==paste(na.s,na.s,sep=na.sep)|is.na(x)|x==na.s,na.s=na.strings,na.sep=allele.sep)
  .data <- apply(.data,2,as.character)
  .data <- replace(.data,.data==na.strings,paste(na.strings,na.strings,sep=allele.sep))
  ## IF THE LENGTH OF THE CHARACTER IS GREATER THAN 2:
  if(all(nchar(.data[.true == F]) > 2))
    {
      .data <-lapply(1:dim(.data)[2],function(i,data,allele.sep) do.call("rbind",strsplit(unlist(data[,i]),split=allele.sep)),data=.data,allele.sep=allele.sep)
      .data <- do.call("cbind",.data)
    }
  else
    {
      if(all(nchar(.data[.true == F]) == 2))
        {
          .data <- replace(.data,.data == paste(na.strings,allele.sep,na.strings,sep=""),paste(na.strings,na.strings,sep=""))
          
          if(T){ ## SPLIT ALL CHARACTER STRINGS IN HALF, PUT SIDE BY SIDE IN SEPARATE COLUMNS
		.dim.tmp <- dim(.data)
		.nchar <- nchar(.data)
		.data <- rbind(substring(.data, first = 1, last = .nchar/2), substring(.data, first = .nchar/2 + 1, last = .nchar))
		dim(.data) <- c(.dim.tmp[1], 2 * .dim.tmp[2])
	}
          else{# GAMMEL FRA HGB
          .mat <- nchar(.data)/2
          .data <- lapply(1:dim(.data)[2],function(i,data,.mat) lapply(1:dim(.mat)[1],function(i,j,data,.mat) paste(substring(data[j,i],1,.mat[j,i]),substring(data[j,i],1+.mat[j,i],.mat[j,i]+.mat[j,i]),sep="."),data = data,.mat=.mat,i=i),data=.data,.mat=.mat)
          .data <- do.call("cbind",.data)
          .data <- lapply(1:dim(.data)[2],function(i,data) do.call("rbind",strsplit(unlist(data[,i]),split="\\.")),data=.data)
          .data <- do.call("cbind",.data)
	}# END GAMMEL	
        }
    }
  ## FIND THE NAMES OF THE UNIQUE ALLELES IN EACH LOCUS:
  .data <- apply(.data,2,as.character)
  ## IF ARGUMENT FAMILY IS EQUAL TO "MFC" THEN THERE IS 6 COLUMNS FOR EACH LOCI, ELSE IF FAMILY IS
  ## EQUAL TO "C" THEN THER IS 2 COLUMNS FOR EACH LOCI
  if(family == "mfc")
    {
      .t <- 6
    }
  else
    {
      if(family != "c") stop(" Problem with family argument!")
      .t <- 2
    }
  ## FIND THE NUMBER OF LOCUS
  .nloci <- (dim(.data)[2])/.t
  if(is.numeric(markers))
    {
      ## PICK OUT THE MARKERS SPECIFIED IN MARKERS ARGUMENT 
      .data <- do.call("cbind",lapply(1:length(markers),function(i,markers,.data,.t) .data[,(1+(markers[i]-1)*.t):(markers[i]*.t)],markers=markers,.data=.data,.t=.t))
      .nloci <- length(markers)
    }
  row.names(.data) <- 1:(dim(.data))[1]
  ## REMOVE MISSING ROWS :
  if(!use.missing)
    {
      ## CAN'T USE NA.EXCLUDE DUE TO DATA ON CHARACTER FORMAT
      if(variables > 0) .xdata <- .xdata[apply(.data,1,function(x,na.strings) all(x != na.strings),na.strings=na.strings),,drop = F]
      .data <- .data[apply(.data,1,function(x,na.strings) all(x != na.strings),na.strings=na.strings),]
      .na.message <- dim(.full.data)[1]-dim(.data)[1]
    }
  else
    {
    	.tmp <- ((.data == "NA") + 0) %*% rep(1,dim(.data)[2]) == 0 # Q & D
      .na.message <- dim(.full.data)[1]-sum(.tmp)
    }
  .row.names <- row.names(.data)
  ## FIND HOW MANY ROWS REMOVED DUE TO MISSING VALUES:
  .omitted <- seq(nrow(.full.data))[is.na(match(row.names(.full.data),row.names(.data)))]
  if(length(.omitted) == 0) .omitted <- NULL
#### VET IKKE HVORFOR DETTE VAR MED:
###  ## CONUNT NUMBER OF LOCI and ALLELES in indata. COUNT DIFFERENT ALLELES IN EACH LOCUS:
###  if(nchar(.nloci) > 1)
###    stop("Problem with number of entries of alleles in markers")
  ## CONSTRUCT A NAMEVECTOR FOR DATA :
  .namevector <- NULL
  if(is.numeric(markers))
    {
      if(family == "mfc")
        {
          .namevector <- c(sapply(1:.nloci,function(i,markers) c(paste("l",markers[i],".m1",sep=""),paste("l",markers[i],".m2",sep=""),paste("l",markers[i],".f1",sep=""),paste("l",markers[i],".f2",sep=""),paste("l",markers[i],".c1",sep=""),paste("l",markers[i],".c2",sep="")),markers=markers))
        }
      else
        {
          .namevector <- c(sapply(1:.nloci,function(i,markers) c(paste("l",markers[i],".c1",sep=""),paste("l",markers[i],".c2",sep="")),markers=markers))
        }
    }
  else
    {
      if(family == "mfc")
        {
          .namevector <- c(sapply(1:.nloci,function(i) c(paste("l",i,".m1",sep=""),paste("l",i,".m2",sep=""),paste("l",i,".f1",sep=""),paste("l",i,".f2",sep=""),paste("l",i,".c1",sep=""),paste("l",i,".c2",sep=""))))
        }
      else
        {
          .namevector <- c(sapply(1:.nloci,function(i) c(paste("l",i,".c1",sep=""),paste("l",i,".c2",sep=""))))
        }
    }
  ## IF THERE ARE ANY VARIABLES IF DATA, THEY ARE ADDED :
  if(variables > 0)
    {
      .namevector <- c(.xnamevec,.namevector)
      .data <- cbind(apply(.xdata[,,drop=F],2,as.character),.data)
    }
  if(use.missing & F)
    {
      non.na.rows <- row.names(na.exclude(.data))
      .match <- match(row.names(.data),non.na.rows)
      na.rows <- seq(nrow(.data))[is.na(.match)]
      .data <- rbind(.data[non.na.rows,],.data[na.rows,])
      row.names(.data) <- c(non.na.rows,na.rows)
    }
  ## MAKE SURE THAT DATA ARE CHARACTER
  .data <- apply(.data,2,as.character)
  dimnames(.data) <- list(NULL,.namevector)
  ## RETURNS NUMBER OF ROWS DROPPED:
  attr(.data,"rows.with.na")  <- .na.message
  ## RETURNS WHICH ROWS DROPPED:
  attr(.data,"rows.dropped")  <- as.numeric(.omitted)
  if(is.numeric(markers)) attr(.data,"markers") <- markers
  return(.data)
}
 
