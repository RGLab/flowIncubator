#' function that runs the dimension-reduction algorithm tSNE (t-Distributed Stochastic Neighbor Embedding, Van der Maaten's Barnes-Hut implementation, R pkg 'Rtsne') on a gatingSet
#' Will sample the minimal number of cells available in all samples to generate balanced cell counts
#' 
#' IMPORTANT: Requires a valid gatingSet with cytokine gates downstream of a parent gate
#' Also expects that pData(gs) contains at least columns: 'name', 'ptid'
#' 
#' @param gs a GatingSet object, properly gated data with annotation in its pData
#' @param parentGate a \code{string} describing the gate upstream of the cytokine gates (eg. "CD4", "cd8+", etc...)
#' @param cytokine a \code{vector} of \code{strings} describing the cytokine gates immediately downstream of parentGate, eg: "IL2", "IFNg"
#' @param otherMarkers the remaining markers of the data
#' @param markerMap named list of marker names to gate names, eg.  list("CD4/IL2" = "IL2","CD4/IFNg" = "IFNg")
#' @param groupBy columns of the \code{gatingSet}
#' @param seed a seed since tSNE is random
#' @param theta parameter to be passed to the \code{Rtsne} function
#' @param ... other parameters to be passed to the \code{Rtsne} function
#' @return a \code{matrix} of X and Y coordinates
#' 
require(flowWorkspace)
require(flowIncubator)
require(data.table)
require(plyr)
require(Rtsne)

runTSNE <- function(gs, parentGate, cytokines, otherMarkers, markerMap, groupBy, seed=999, theta=0.9, ...){
  set.seed(seed) 
  
  #get pData
  pd <- as.data.table(pData(gs))
  
  #get cell counts of the parent gate
  cat("getting total cell counts from parent gate", parentGate, "\n")
  parent_count <- unlist(lapply(gs, function(gh)getTotal(gh, parentGate)))
  parent_count = ldply(parent_count)
  
  # col <- eval(paste(parentGate,"count", sep="_"))
  
  setnames(parent_count,c("name",parentGate))
  
  
  pd <- merge(pd, parent_count,by="name")
  
  
  nTcells <- min(pd[, sum(get(parentGate)), by = groupBy][, V1])
  cat("all samples have at least", nTcells, "cells. This will be the sample size...  ")
  # NOTE: this will fail if one or more sample has no cells !!!
  
  pd[, {
    
    #subsample entire group together
    totalEvents <- sum(get(parentGate))
    gInd <- 1:totalEvents
    gInd <- sample.int(totalEvents, size = nTcells)
    #convert it to logical vector
    gInd.logical <- rep(F, totalEvents)
    gInd.logical[gInd] <- T
    
    #split it into sample level
    sn.factor <- unlist(sapply(name, function(sn)rep(sn, .SD[name == sn, get(parentGate)])))
    ind.vec <- split(gInd.logical, sn.factor)
    #       
#     #reset indices for each sample
    for(sn in name)
    {
      thisInd <- ind.vec[[sn]]
      gh <- gs[[sn]]
      updateIndices(gh, parentGate, thisInd)
    }
    
  }
  , by = groupBy]
  cat("resampling complete ! recomputing... ")
  
  nodes <- getChildren(gs[[1]], parentGate, path = 2)
  
  #' recompute down-stream gates due to the resetting of parent gate indices
  for(node in nodes)
    recompute(gs, node)
  
  res <- getSingleCellExpression(gs
                                 , nodes
                                 , other.markers = otherMarkers
                                 , map = markerMap
                                 , threshold = FALSE
  )
  
  res_mask <- getSingleCellExpression(gs
                                      , nodes
                                      ,  map = markerMap
                                      , threshold = T
  )
  
  #' collapse the single Cell expression data matrices across samples
  res_collapse <- ldply(names(res), function(sn){
    message(".", appendLF = F)
    mat <- as.data.frame(res[[sn]])
    
    if(nrow(mat) > 0){
      #compute degree of functionality      
      mat_mask <- res_mask[[sn]]
      mat_mask[mat_mask > 0] <- 1
      mat[, "degree"] <- rowSums(mat_mask)
      #compute the boolean expression of each event
      mat[, "poly"] <- Reduce(paste0, as.list(as.data.frame(mat_mask)))
      pd <- pData(gs[[sn]])
      rownames(pd) <- NULL
      cbind(mat, pd)  
    }else
      NULL
  })
  
  cat("\n cytokines:", cytokines)
  cat("\n other markers:", otherMarkers)
  
  included_markers <- c(cytokines, otherMarkers)
  input_mat <- as.matrix(res_collapse[,included_markers])
  
  cat("input is ready, starting tSNE run... \n")
  system.time(tsne_out <- Rtsne(input_mat, check_duplicates = FALSE, ...))
  
  dat <- tsne_out$Y
  colnames(dat) <- c("x", "y")
  dat <- cbind(dat, res_collapse)
  dat <- data.table(dat)
  
  cat("DONE !")
  return(dat)
  
}