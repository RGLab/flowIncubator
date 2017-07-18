#' impute the gate (flagged as failure by external algorithm) with refGate from nearest neighbour sample
#'
#' @param gs a \code{GatingSet}
#' @param node a \code{character} or \code{numeric} specifing node index
#' @param failed a \code{character} or \code{numeric} specifing samples that fails the gating QA
#' @param passed a \code{character} or \code{numeric} specifing samples that passes QA to be served as references
#'                                                        By default, all samples other than failed will be used.
#'                                                        but sometime it is helpful to narrow it down to a few of really good samples. 
#' @param ... other arguments passed to \code{.nearestSample}
#'            n an \code{integer} passed to \code{density} call. default is 512.
#'            bandwidth, gridsize passed to \code{bkde2D} call. default is c(5, 5) and c(100, 100).
#'            method density similarity calculation method. Either "ks.test" or "em". For 2d gate, "em" is the only choice.
#'            mc.cores passed to mclapply for parallel computing. Default is 1.
#' @examples 
#' \dontrun{
#' library(flowWorkspace)
#' gs <- load_gs("gs-tcell/")
#' 
#' #Finding reference sample for: 1349_3_Tcell_A06.fcs
#' res <- nearestSamples(gs, node = "CD3", failed = "1349_3_Tcell_A06.fcs")
#' # res is a named vector stores the mapping between bad(failed) and good(reference) samples
#' 
#' #fix the gate for the bad samples
#' regateNearestSamples(gs, res, "CD3")
#' 
#' #run it on 2d gate and customize some parameters
#' #and enable parallel computing to speed it up
#' res <- nearestSamples(gs, node = "CD4", failed = "1349_3_Tcell_A06.fcs", gridsize = c(70, 70), mc.cores = 4)
#' #fix the 2d gate for the bad samples
#' regateNearestSamples(gs, res, "CD4")

#' }
#' @export
nearestSamples <- function(gs, node, failed, passed = NULL, ...){
  #get samples that do not fail the QA check
  failed <- as.vector(failed)
  samples <- sampleNames(gs)
  failedInd <- match(failed,samples)
  if(is.null(passed))
    passed <- samples[-c(failedInd)]
  
  sapply(failed,function(thisTarget){
    message("Finding reference sample for: ",thisTarget)
    .nearestSample(gs, node = node, target = thisTarget, source = passed, ...) 
  })
  
}

#' find the nearest neighbour sample that shares the most similar density profile of the specifed gate   
#' the similarity is defined by the Earth-Mover's distance 
#'
#' @param gs a \code{GatingSet}
#' @param node a \code{character} or \code{numeric} specifing node index
#' @param n an \code{integer} passed to \code{density} call
#' @param ... other arguments passed to \code{density} call
#' @importFrom emdist emd emd2d
#' @importFrom KernSmooth bkde2D
#' @importFrom parallel mclapply
.nearestSample <- function(gs, node, target, source, n = 512, bandwidth = c(5, 5), gridsize = c(100, 100), method = c("ks.test","em"), mc.cores = 1,  ...){
  method <- match.arg(method)
  
  thisGh <- gs[[target]]
  
  #get data 
  parentNode <- getParent(thisGh, node)
  thisGate <- getGate(thisGh, node)
  params <- parameters(thisGate)
  nDim <- length(params)
  parentData <- getData(gs,parentNode)[,params]
  #get data
  tData <- parentData[[target]]
  
  #exclude marginal events below zero
  #    expression1 <- paste0("`",params,"`>0")
  #    ef <- char2ExpressionFilter(expression1)
  #    tData <- Subset(tData,ef)
  
  tExpr <- exprs(tData)
  
  if(nDim > 2)
    stop("Imputing Gate on more than 2 dimensions is not supported yet!")
  
  if(method == "em"){
    #get density of failed sample
    if(nDim == 1)
    {
      tDen <- density(tExpr, n =n ,...)
      tMat <- matrix(c(tDen$y,tDen$x),ncol = 2)  
    }else
    {
      tDen <- bkde2D(tExpr, bandwidth = bandwidth,  gridsize = gridsize, ...)
      tMat <- tDen$fhat
      # contour(tDen$x1, tDen$x2, tDen$fhat)
    }
    
  }
  #TODO:customize mc.cores
  #cal dist from each sample
  thisCall <- quote(lapply(source,function(thisSample){
    #get 1d density of target sample
    thisData <- parentData[[thisSample]]
    #apply the same filter
    #                              thisData <- Subset(thisData,ef)
    
    thisExpr <- exprs(thisData)
    #cal the dist
    #                              browser()
    if(method == "ks.test"){
      thisDist <- sapply(seq_len(nDim), function(i){
        ksstat(tExpr[,i],thisExpr[,i])
      })
      thisDist <- sum(thisDist)
    }else if(method == "em"){
      #EM
      if(nDim == 1){
        thisDen <- density(thisExpr, n =n ,...)
        thisMat <- matrix(c(thisDen$y,thisDen$x),ncol = 2)
        thisDist <- emd(tMat,thisMat)    
      }else
      {
        thisDen <- bkde2D(thisExpr, bandwidth = bandwidth,  gridsize = gridsize, ...)
        # contour(thisDen$x1, thisDen$x2, thisDen$fhat)
        thisMat <- thisDen$fhat
        thisDist <- emd2d(tMat,thisMat)    
      }
      
    }
      
      thisDist
  }))
  
  if(mc.cores > 1)
  {
    thisCall[[1]] <- quote(mclapply)
    thisCall[["mc.cores"]] <- mc.cores
  }
    
  distVec <- eval(thisCall)
    #    browser()
    ##visualize the distance vs density
    #    dat <- parentData[c(target,source[-3])]
    #    dat <- Subset(dat,ef)
    #    pData(dat)$dist <- c(0,as.numeric(format(unlist(distVec[-3]),digits=2)))
    #    densityplot(as.formula(paste0("as.factor(dist)~`",params,"`")),dat
    #            ,darg=list(bw="nrd0",n=n)
    #              , main = paste(method, "distance")
    #            )
    #pick the closet one    
    source[which.min(distVec)]
    
    
  
    
  
}

#'  perform the actual regating of identified bad samples with gates from good samples
#'  using the named vector from flowIncubator:::.nearestSamples()
#'
#' @param gs a \code{GatingSet}
#' @param samples a named \code{character} vector of sample names in the supplied GatingSet
#' @param gate_name a \code{character} specifying the gate to operate on (eg. "CD3")
#' @export
regateNearestSamples <- function (gs, samples, gate_name){
  for (i in 1:length(samples)){
    bad <- names(samples)[i] 
    good <- getGate(gs[samples[i]], gate_name)    
    cat("regating", gate_name, "replacing", bad, "with", samples[i], "...")
    names(good) <- bad
    setGate(gs[bad], gate_name, good)
    recompute(gs[bad], gate_name)
  }  
}
