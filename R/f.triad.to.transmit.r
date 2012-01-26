f.triad.to.transmit <-
  function(indata,outdata,sep = " ",allele.sep=";",na.strings="NA")
{
  ## INDATA IS AN ASCII FILE WHERE DATA IS STORED. SEP IS THE 
  ## SEPARATOR BETWEEN COLUMNS IN INDATA.
  ## THIS FUNCTION CONVERTS DATA LISTED AS TRIADE FORMAT TO TRANSMIT FORMAT.
  ## THE ARGUMENT SEP INDICATE THE BETWEEN COLUMN SEPARATOR WHILE THE ALLELE.SEP
  ## INDICATE THE WITHIN COLUMN SEPARATOR:
  indata <- as.character(indata)
  if(allele.sep==".") allele.sep <- "\\."
  ## IF THE SEPARATOR IS SOME SORT OF WHITE SPACE, THE SEP ARGUMENT IN COUNT.FIELDS AND SCAN ARE LEFT OUT:
   if(na.strings==" " || na.strings=="")
     {
       ## COUNT FIELDS, SET UP NAMES:
       .count <- count.fields(file = indata,sep=sep)
       .nvar <- .count[1]
       if(any(.count != .nvar))
         stop("Number of variables differ from row to row")
       .data <- read.table(file = indata,sep=sep,na.strings=na.strings, stringsAsFactors = F)
     }
   else
     {
       ## COUNT FIELDS, SET UP NAMES:
       .count <- count.fields(file = indata)
       .nvar <- .count[1]
       if(any(.count != .nvar))
         stop("Number of variables differ from row to row")
       .data <- read.table(file = indata,na.strings=na.strings, stringsAsFactors = F)
     }
  .data <- lapply(1:length(.data),function(i,.data) as.character(.data[[i]]),.data=.data)
  .data <- lapply(1:length(.data),function(i,.data,allele.sep) replace(.data[[i]],.data[[i]]=="NA",paste("NA","NA",sep=allele.sep)),.data=.data,allele.sep=allele.sep)
  .data <- do.call("cbind",.data)
  ## IF THE DOES NOT EXIST ANY NA.STRING VALUES IN .DATA:
  if(all(nchar(.data) > 2))
    {
      .data <-lapply(1:dim(.data)[2],function(i,data,allele.sep) do.call("rbind",strsplit(unlist(data[,i]),split=allele.sep)),data=.data,allele.sep=allele.sep)
      .data <- do.call("cbind",.data)
    }
  if(all(nchar(.data) == 2))
    {
      .data <- lapply(1:dim(.data)[2],function(i,data) paste(substring(data[,i],1,1),substring(data[,i],2,2),sep="."),data=.data)		
      
      .data <- do.call("cbind",.data)
      .data <- lapply(1:dim(.data)[2],function(i,data) do.call("rbind",strsplit(unlist(data[,i]),split="\\.")),data=.data)			
      .data <- do.call("cbind",.data)	
    }
  .nloci <- (dim(.data)[2])/6
  .nallele <- NULL
  .unique <- lapply(1:.nloci,function(i,data) names(table(data[,(1+(i-1)*6):(6*i)])),data=.data)
  .nallele <- as.vector(unlist(lapply(1:.nloci,function(i,x) length(x[[i]]),x=.unique)))
  .sum <- lapply(1:.nloci,function(i,data) table(data[,(1+(i-1)*6):(6*i)],exclude=NULL),data=.data)
  ## PRINTS NUMBER OF LOCI AND THE RESPECTIVE NUMBER OF ALLELES:
  cat("Number of marker(s):",.nloci,"\n")
  unlist(sapply(1:.nloci,function(i, uni,summ,nallele) cat("Marker",i,"has",nallele[i],"different alleles","\n",paste("Allele",as.vector(unlist(uni[[i]])),"has frequency",as.vector(unlist(summ[[i]])),"\n"),"\n"),uni=.unique,summ=.sum,nallele=.nallele))	
  .namevector <- NULL
  .namevector <- c(sapply(1:.nloci,function(i) c(paste("l",i,".m1",sep=""),paste("l",i,".m2",sep=""),paste("l",i,".f1",sep=""),paste("l",i,".f2",sep=""),paste("l",i,".c1",sep=""),paste("l",i,".c2",sep=""))))	
  dimnames(.data) <- list(NULL,.namevector)
  ## CONVERTING DATA TO TRANSMIT FORMAT:
  .nrow <- dim(.data)[1]
  .names <- c("family","id","father","mother","sex","affected",paste("marker", 1.:.nloci, sep = ""))
  .family <- sort(rep(seq(1:.nrow),3))
  .id <- rep(1:3,.nrow)
  .father <- rep(c(".",".","1"),.nrow)
  .mother <- rep(c(".",".","2"),.nrow)
  .sex <- rep(c("1","2","0"),.nrow)
  .affected <- rep(c("0","0","0"),.nrow)
  .namefather <- lapply(1:.nloci,function(i) c(paste("l",i,".f1",sep=""),paste("l",i,".f2",sep="")))
  .namemother <- lapply(1:.nloci,function(i) c(paste("l",i,".m1",sep=""),paste("l",i,".m2",sep="")))
  .namechild <- lapply(1:.nloci,function(i) c(paste("l",i,".c1",sep=""),paste("l",i,".c2",sep="")))	
  .marker.f <- lapply(1:.nloci,function(i,data,namefather) paste(data[,namefather[[i]][1]],"/",data[,namefather[[i]][2]],sep=""),data=.data,namefather=.namefather)
  .marker.m <- lapply(1:.nloci,function(i,data,namemother) paste(data[,namemother[[i]][1]],"/",data[,namemother[[i]][2]],sep=""),data=.data,namemother=.namemother)
  .marker.c <- lapply(1:.nloci,function(i,data,namechild) paste(data[,namechild[[i]][1]],"/",data[,namechild[[i]][2]],sep=""),data=.data,namechild=.namechild)	
  .marker <- lapply(1:.nloci,function(i,marker.f,marker.m,marker.c) as.vector(do.call("rbind",list(marker.f[[i]],marker.m[[i]],marker.c[[i]]))),marker.f=.marker.f,marker.m=.marker.m,marker.c=.marker.c)
  .marker <- do.call("cbind",.marker)
  .transmit <- do.call("cbind",list(.family,.id,.father,.mother,.sex,.affected,.marker))
  dimnames(.transmit) <- list(NULL,.names)
  ## ADD EXTRA LINE AT THE TOP WHICH INCLUDE THE NUMER OF MARKER IN THE DATASET:
  .transmit <- rbind(c(.nloci,rep(" ",dim(.transmit)[2]-1)),.transmit)
  ## WRITE RESULT TO OUTDATA:
  write.table(.transmit, file = outdata, sep = " ", append = F, quote = FALSE,row.names = FALSE,col.names = FALSE,na="NA")
}
