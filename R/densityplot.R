#' plot the population as densityplot for the channels associated with the gate
#' 
#' It is used to check the 1d density of the parent for the purpose of choosing appropriate 
#' 1d gating algorithm or fine tune the gating parameters.
#'  
#' @param x GatingHierarchy
#' @param data node/population name
setMethod("densityplot", signature = c("GatingHierarchy", "character")
          , definition = function(x, data, ...){
  gh <- x
  node <- data
  parent <- getParent(gh, node)
  parent <- getData(gh, parent)
  gate <- getGate(gh, node)
  chnls <- parameters(gate)
  densityplot(~., parent, channels = chnls, ...)
})