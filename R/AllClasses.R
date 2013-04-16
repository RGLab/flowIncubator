
setClass("GatingSetList"
        ,representation=representation(
            data = "list"
            ,samples="character" #this determine the order of samples exposed to user
          )
        )   
## Constructor
GatingSetList <- function(x)
{
  names(x)<-NULL#strip names from the list because rbind2 doesn't like it
  flowCore:::checkClass(x, "list")
  
  x <- new("GatingSetList", data = x, samples = unlist(lapply(x,getSamples)))
  return(x)
}