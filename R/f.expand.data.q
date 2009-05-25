f.expand.data <- function(data)
  {
    .namevector <- NULL
    .data <- data
    .alleles <- attr(.data,"alleles")
    .nloci <- length(.alleles)
    ## FIND IF FAMILY ARGUMENT FROM F.READ.DATA IS "MFC" OR "C"
    if(any(regexpr("\\.m",dimnames(data)[[2]]) > 0))
      {
        .family <- "mfc"
        .t <- 4
      }
    else
      {
        .family <- "c"
        .t <- 2
      }
    if(!is.null(attr(.data,"markers")))
      {
        .markers <- attr(.data,"markers")
        .nloci <- length(.markers)
        .allele.names <- lapply(1:.nloci,function(i,.alleles,.markers) match(names(.alleles[[i]]),sort(unique.default(names(.alleles[[i]])), na.last = T)),.alleles=.alleles,.markers=.markers)
      }
    else
      {
        .markers <- 1:.nloci
        .allele.names <- lapply(1:.nloci,function(i,.alleles) match(names(.alleles[[i]]),sort(unique.default(names(.alleles[[i]])), na.last = T)),.alleles=.alleles)
      }
    .temp.allelelist <- rep(lapply(1:.nloci, function(i, .allele.names) as.numeric(.allele.names[[i]]), .allele.names=.allele.names), rep(.t, .nloci))
    if(.family == "mfc")
      {
        .namevector <- c(sapply(1:.nloci,function(i,.markers) c(paste("l",.markers[i],".m1",sep=""),paste("l",.markers[i],".m2",sep=""),paste("l",.markers[i],".f1",sep=""),paste("l",.markers[i],".f2",sep="")),.markers=.markers))
      }
    else
      {
        .namevector <- c(sapply(1:.nloci,function(i,.markers) c(paste("l",.markers[i],".c1",sep=""),paste("l",.markers[i],".c2",sep="")),.markers=.markers))
      }
    if(!is.null(attr(.data,"variables")))
      {
        .variables <- attr(.data,"variables")
        .variables <- lapply(1:length(.variables),function(i,.variables) match(names(.variables[[i]]),sort(unique.default(names(.variables[[i]])), na.last = T)),.variables=.variables)
        .temp.allelelist <- c(.variables,.temp.allelelist)
        .xnamevec <- c(sapply(1:length(attributes(.data)$variables),function(i) c(paste("x",i,sep=""))))
      }
    if(!is.null(attr(.data,"variables"))) .namevector <- c(.xnamevec,.namevector)
    ## EXPAND ENTIRE GRID
    .grid <- do.call("expand.grid", .temp.allelelist)
    ## CREATE ARRAY 
    .arr <- array(data = 0,dim = sapply(.temp.allelelist,length))
    .data <- .data[,c(.namevector,"freq"),drop=F]
    .index <- .data[,1:(dim(.data)[2]-1),drop=F]
    .arr[.index] <- .data[,"freq"]
    .grid$freq <- as.numeric(.arr)
    dimnames(.grid)[[2]] <- c(.namevector,"freq")
    ## RETURNS ALLELENAMES AS AN ATTRIBUTE TO .DATA:
    attr(.grid,"alleles") <- .alleles
    if(!is.null(attr(.data,"variables"))) attr(.grid,"variables") <- attr(.data,"variables")
    return(.grid)
  }

