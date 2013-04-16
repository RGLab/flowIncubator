
setClass("GatingSetList"
        ,representation=representation(
            data = "list")
        )   
## Constructor
GatingSetList <- function(x)
{
  names(x)<-NULL#strip names from the list because rbind2 doesn't like it
  flowCore:::checkClass(x, "list")
  x <- new("GatingSetList", data=x)
  return(x)
}