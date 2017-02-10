#' run tSNE from (R pkg 'Rtsne') on a gatingSet
#' Will sample the minimal number of cells available in all samples to generate balanced cell counts
#' 
#' IMPORTANT: Requires a valid gatingSet with cytokine gates downstream of a parent gate
#' Also expects that pData(gs) contains at least columns: 'name', 'ptid' so we can identify cells later
#' 
#' @param gs a GatingSet object, properly gated data with annotation in its pData
#' @param parentGate a \code{string} describing the gate upstream of the cytokine gates (eg. "CD4", "cd8+", etc...)
#' @param cytokine a \code{vector} of \code{strings} describing the cytokine gates immediately downstream of parentGate, eg: "IL2", "IFNg"
#' @param otherMarkers the remaining markers of the data
#' @param markerMap named list of marker names to gate names, eg. list("CD4/IL2" = "IL2","CD4/IFNg" = "IFNg")
#' @param swap boolean for whether marker and gate names (from markerMap above) should be swapped. Passed onto getSingleCellExpression()
#' @param groupBy columns of the \code{gatingSet}'s phenoData, same number of cells will be sampled from each group
#' @param degreeFilter keep cells of this degree and higher, useful when tSNE takes too long to run
#' @param seed since tSNE is random, need a random seed so we can reproduce results
#' @param theta parameter to be passed to the \code{Rtsne} function
#' @param ... other parameters to be passed to the \code{Rtsne} function
#' @return a \code{matrix} of X and Y coordinates
#' @import flowWorkspace
#' @import data.table
#' @import plyr
#' @import Rtsne
runTSNE <- function (gs, parentGate, cytokines, otherMarkers, markerMap, swap = FALSE,
                     groupBy, degreeFilter = 0, seed = 999, theta = 0.9, ...) {
  
  if (is.null(markerMap)) stop ("required markerMap is missing ! STOPPING....")
    
  set.seed(seed)
  pd <- as.data.table(pData(gs))
  meta_cols <- colnames(pd)
  meta_cols <- c(meta_cols, "degree", "poly")
  cat("getting total cell counts from parent gate", parentGate, 
      "\n")
  parent_count <- unlist(lapply(gs, function(gh) getTotal(gh, 
                                                          parentGate)))
  parent_count = ldply(parent_count)
  setnames(parent_count, c("name", parentGate))
  pd <- merge(pd, parent_count, by = "name")
  nTcells <- min(pd[, sum(get(parentGate)), by = groupBy][, 
                                                          V1])
  cat("after grouping by '", groupBy, "', all groups will have at least", 
      nTcells, "cells.\n")
  pd[, {
    totalEvents <- sum(get(parentGate))
    gInd <- 1:totalEvents
    gInd <- sample.int(totalEvents, size = nTcells)
    gInd.logical <- rep(F, totalEvents)
    gInd.logical[gInd] <- T
    sn.factor <- unlist(sapply(name, function(sn) rep(sn, 
                                                      .SD[name == sn, get(parentGate)])))
    ind.vec <- split(gInd.logical, sn.factor)
    for (sn in name) {
      thisInd <- ind.vec[[sn]]
      gh <- gs[[sn]]
      updateIndices(gh, parentGate, thisInd)
    }
  }, by = groupBy] 
  cat("subsampling complete ! recomputing... \n")
  nodes <- getChildren(gs[[1]], parentGate, path = 2)
  for (node in nodes) recompute(gs, node)
  cat("generating event masks \n")
  
  res <- getSingleCellExpression(gs, nodes, other.markers = otherMarkers, 
                                 map = markerMap, threshold = FALSE, swap=swap)
  
  res_mask <- getSingleCellExpression(gs, nodes, map = markerMap, 
                                      threshold = T, swap=swap)
  
  res_collapse <- ldply(names(res), function(sn) {
    message(".", appendLF = F)
    # message(sn)
    mat <- as.data.frame(res[[sn]])
    if (nrow(mat) > 0) {
      mat_mask <- res_mask[[sn]]
      mat_mask[mat_mask > 0] <- 1
      mat[, "degree"] <- rowSums(mat_mask)
      mat[, "poly"] <- Reduce(paste0, as.list(as.data.frame(mat_mask)))
      pd <- pData(gs[[sn]])
      rownames(pd) <- NULL
      cbind(mat, pd)
    } 
    else NULL
  })
  
  
  cat("\n cytokines:", cytokines)
  cat("\n other markers:", otherMarkers)
  cat("\n input matrix has", nrow(res_collapse), "rows...")
  res_collapse <- subset(res_collapse, degree > degreeFilter)
  cat("\n input matrix has", nrow(res_collapse), "rows after filtering for cells of degree >", 
      degreeFilter)
  input_mat <- as.matrix(res_collapse[, !names(res_collapse) %in% 
                                        meta_cols])
  cat("\n starting tSNE run at ", date(), "\n")
  system.time(tsne_out <- Rtsne(input_mat, check_duplicates = FALSE, 
                                ...))
  dat <- tsne_out$Y
  colnames(dat) <- c("x", "y")
  dat <- cbind(dat, res_collapse)
  dat <- data.table(dat)
  cat("completed tSNE run at", date(), "!\n")
  return(dat)
  
}