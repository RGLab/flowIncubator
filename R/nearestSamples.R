#' impute the gate (flagged as failure by external algorithm) with refGate from nearest neighbour sample
#'
#' @param gs a \code{GatingSet}
#' @param node a \code{character} or \code{numeric} specifing node index
#' @param failed a \code{character} or \code{numeric} specifing samples that fails the gating QA
#' @param passed a \code{character} or \code{numeric} specifing samples that passes QA to be served as references
#'                                                        By default, all samples other than failed will be used.
#'                                                        but sometime it is helpful to narrow it down to a few of really good samples. 
#' @param ... other arguments passed to \code{.nearestSample}
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
#' }
#' @export
nearestSamples <- function(gs, node, failed, passed = NULL, ...){
  #get samples that do not fail the QA check
  #  browser()
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
.nearestSample <- function(gs, node, target, source, n = 512, method = c("ks.test","em"), ...){
  method <- match.arg(method)
  
  thisGh <- gs[[target]]
  
  #get data 
  parentNode <- getParent(thisGh, node)
  thisGate <- getGate(thisGh, node)
  params <- parameters(thisGate)
  parentData <- getData(gs,parentNode)[,params]
  
  if(length(params) == 1){
    
    #get data
    tData <- parentData[[target]]
    
    #exclude marginal events below zero
    #    expression1 <- paste0("`",params,"`>0")
    #    ef <- char2ExpressionFilter(expression1)
    #    tData <- Subset(tData,ef)
    
    tExpr <- exprs(tData)
    
    if(method == "em"){
      #get 1d density of failed sample
      tDen <- density(tExpr, n =n ,...)
      tMat <- matrix(c(tDen$y,tDen$x),ncol = 2)
    }
    #TODO:customize mc.cores
    #cal dist from each sample
    distVec <- lapply(source,function(thisSample){
      #get 1d density of target sample
      thisData <- parentData[[thisSample]]
      #apply the same filter
      #                              thisData <- Subset(thisData,ef)
      
      thisExpr <- exprs(thisData)
      #cal the dist
      #                              browser()
      if(method == "ks.test"){
        thisDist <- ks.test(tExpr,thisExpr)$statistic
      }else if(method == "em"){
        #EM
        
        thisDen <- density(thisExpr, n =n ,...)
        thisMat <- matrix(c(thisDen$y,thisDen$x),ncol = 2)
        thisDist <- emd(tMat,thisMat)  
      }
      ks.test <-
        
        
        thisDist
    })
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
    
    
  }else if(length(params) == 2){
    stop("Imputing 2d Gate not supported yet!")
  }
  
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
