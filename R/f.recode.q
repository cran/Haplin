
f.recode <- function(data.new,data.old,ind)
  {
    for(i in 1:length(ind))
      {
        data.new[data.old == ind[i]] <- i
      }
    return(data.new)
  }
