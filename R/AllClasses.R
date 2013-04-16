
setClass("GatingSetList"
        ,representation=representation(
            data = "list"
            ,samples="character" #this determine the order of samples exposed to user
          )
        )   
## Constructor
GatingSetList <- function(x,samples = NULL)
{
  names(x)<-NULL#strip names from the list because rbind2 doesn't like it
  flowCore:::checkClass(x, "list")
  if(is.null(samples)){
    samples <- unlist(lapply(x,getSamples))
  }
  x <- new("GatingSetList", data = x, samples = samples)
  return(x)
}