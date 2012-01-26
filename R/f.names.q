f.names <-  function(i,x,y)
  {
    x <- as.vector(unlist(x[[i]]))
    names(x) <- as.vector(unlist(y[[i]]))
    return(x)
  }
                            
