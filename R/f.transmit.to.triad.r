f.transmit.to.triad <- function(indata,outdata)
{
  ## READS AND CONVERTS DATA FROM FORMAT USED BY transmit TO TRIAD FORMAT.
  ## INDATA IS THE FILE CONTAINING transmit DATA. OUTDATA IS THE 
  ## ASCII FILE OF THE CONVERTED DATA.
  
  ## COUNT FIELDS, SET UP NAMES
  cat("Converting data....\n\n")
  .count <- count.fields(file = indata, sep = " ")
  ## THE FIRST LINE IN TRANSMIT FORMAT CONTAIN THE NUMER OF MARKERS IN INDATA.
  ## OMIT THIS LINE:
  .count <- .count[-1]
  .nlines <- length(.count)
  .nvar <- .count[1.]
  .nmarkers <- .nvar - 6.
  if(.nmarkers <= 0.)
    stop("Not enough variables in file!")
  if(any(.count != .nvar))
    stop("Number of variables differ from row to row")
  .names <- c("family", "id", "father", "mother", "sex", 
              "affected")
  .names <- c(.names, paste("marker", 1.:(.nvar - 6.), sep = ""))
  .fields <- as.list(character(.nvar))
  ## SCAN DATA
  names(.fields) <- .names
  .indata <- scan(file = indata, what = .fields, sep = " ", 
                  multi.line = F, flush = F,skip=1)
  ## CHECK family VARIABLE
  if(any(table(.indata$family, exclude = NULL) != 3.))
    stop("Incorrect family size")
  .nfam <- .nlines/3.
  .famcodes <- unique(.indata$family)
  ## CHECK id VARIABLE	
  if(any(tapply(.indata$id, .indata$family, FUN = function(x)
                length(unique(x))) != 3.)) stop("Incorrect id variable")
  ## CHECK AND RECOMPUTE father AND mother VARIABLES
  .father.pos <- .mother.pos <- logical(.nlines)
  .f.parent <- function(id, parent)
    {
      ## PICK WHAT id CORRESPONDS TO FATHER AND MOTHER 
      .parent <- unique(parent)
      .pos <- is.element(id, .parent)
      if(sum(.pos) != 1.)
        stop("Problem with father/mother codes")
      return(.pos)
    }
  for(i in 1.:.nfam) {
    .ind <- .indata$family == .famcodes[i]
    .father.pos[.ind] <- .f.parent(id = .indata$id[.ind],
                                   parent = .indata$father[.ind])
    .mother.pos[.ind] <- .f.parent(id = .indata$id[.ind],
                                   parent = .indata$mother[.ind])
  }
  if(any(.father.pos & .mother.pos))
    {
      stop("Sorry, father and mother cannot have the same id!")
    }
  .child.pos <- !.father.pos & !.mother.pos
  .matrix <- do.call("cbind",.indata)
  .matrix <- cbind(.matrix,.mother.pos,.father.pos,.child.pos)
  .n <- c(paste("marker", 1.:(.nvar - 6.), sep = ""))
  ## PIC OUT MARKER 1-K FOR MOTHER, FATHER AND CHILD.(TRIADE FORMAT)
  .data <- NULL
  .mother <- NULL
  .father <- NULL
  .child <- NULL
  for(i in 1:length(.n))
    {
      .mother <- cbind(.mother,.matrix[.mother.pos==T,.n[i]])
      .father  <- cbind(.father,.matrix[.father.pos==T,.n[i]])
      .child <- cbind(.child,.matrix[.child.pos==T,.n[i]])
    }
  ## WRITE RESULT TO FILE IN TRIAD FORMAT:
  .mother <-  do.call("rbind",lapply(1:dim(.mother)[1],function(i,.mother) sapply(1:dim(.mother)[2],function(j,.mother,i) paste(unlist(strsplit(.mother[i,j],split="/")),collapse=";"),.mother=.mother,i=i),.mother=.mother))
  .father <- do.call("rbind",lapply(1:dim(.father)[1],function(i,.father) sapply(1:dim(.father)[2],function(j,.father,i) paste(unlist(strsplit(.father[i,j],split="/")),collapse=";"),.father=.father,i=i),.father=.father)) 
  .child <-  do.call("rbind",lapply(1:dim(.child)[1],function(i,.child) sapply(1:dim(.child)[2],function(j,.child,i) paste(unlist(strsplit(.child[i,j],split="/")),collapse=";"),.child=.child,i=i),.child=.child))
  .data <- lapply(1:dim(.mother)[2],function(i,mother,father,child) cbind(mother[,i],father[,i],child[,i]),mother=.mother,father=.father,child=.child)
  .data <- do.call("cbind",.data)
  write.table(.data, file = outdata, sep = " ", append = FALSE, quote = FALSE,row.names = FALSE,col.names = FALSE,na="NA")
}

