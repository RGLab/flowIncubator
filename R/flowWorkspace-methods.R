
swapChannelMarker <- function(gs){
  fs <- getData(gs)
  
  for(sn in sampleNames(fs)){
    fr <- fs@frames[[sn]]
    fs@frames[[sn]] <- swapChannelMarker_flowframe(fr, use.exprs = FALSE)
  }
  
  newColNames <- colnames(fs@frames[[sn]])
  colnames(fs) <- newColNames #assuming the order of colnames between fr and fs were consistent
  flowData(gs) <- fs
  
  fr <- fs@frames[[sn]]
  pd <- pData(parameters(fr))
  pd <- pd[!is.na(pd$desc), 2:1]
  colnames(pd) <- c("old", "new")
#  browser()
  gs <- updateGateParameter(gs, pd)
  flowData(gs) <- fs
  gs
}

#' Preprocesses a Cytotrol flowFrame object
#'
#' Our goal here is to use swap the marker names and the channel names within a
#' \code{flowFrame} object to ensure that the \code{flowFrame} objects across
#' centers can be merged into a single \code{flowSet}.
#'
#'
#' @param fr the \code{flowFrame} object to preprocess
#' @return the updated \code{flowFrame} object containing only the markers of
#' interest
swapChannelMarker_flowframe <- function(fr, use.exprs = TRUE) {
  
  
  fr_rownames <- rownames(pData(parameters(fr)))
  
  # Preprocesses each of the columns in the flow_frame
  for (j in seq_len(length(colnames(fr)))) {
    
    marker_idx <- paste0(fr_rownames[j], "S")
    channel_idx <- paste0(fr_rownames[j], "N")
    
    marker <- description(fr)[[marker_idx]]
    channel <- description(fr)[[channel_idx]]
    
    # In the case the marker name is given, we swap the marker and channel
    # names.
    if (!is.null(marker)) {
      # Converts the marker names to a common name
      marker <- as.vector(marker)
      
      # Updates the channel with the marker
      description(fr)[[channel_idx]] <- marker 
      pData(parameters(fr))[j, "name"] <- marker
      
      # Updates the marker information in the flow_frame with the channel
      description(fr)[[marker_idx]] <- channel
      pData(parameters(fr))[j, "desc"] <- channel
    }
  }
  
  if(use.exprs)
    colnames(exprs(fr)) <- colnames(fr)
  
  # Subset to markers of interest
  fr
}


#' update the gate parameters for a GatingSet
#' 
#' It actually reconstructs a new GatingSet by 
#' copying all the gates with their parameters changed 
#' based on given mapping between the old and new channel names.
#' 
#' 
#' @param gs \code{GatingSet} to work with
#' @param map \code{data.frame} contains the mapping between old and new channel names
#' @return  a new \code{GatingSet} object with the new gate added but share the same flow data with the input 'GatingSet'
updateGateParameter <- function(gs, map){
  
  if(!identical(colnames(map), c("old", "new")))
    stop("'map' must contain two columns: 'old' and 'new'!")
  
  #copy the entire tree structure
  message("cloning tree structure...")
  clone <- gs
  clone@pointer <- .Call("R_CloneGatingSet",gs@pointer,sampleNames(gs))
  #clear the tree
  nodes <- getNodes(gs)
  toRm <- getChildren(gs[[1]], "root")
  for(thisRm in toRm)Rm(thisRm, clone)
  
  nodesToadd <- nodes[-1]
  
  for(node in nodesToadd)
  {
        #copy the other nodes to its parent
        thisParent <- getParent(gs, node)
        popName <- basename(node)
        
        for(sn in sampleNames(gs))
        {
              gh <- gs[[sn]]
              gate <- getGate(gh, node)
              
              if(!flowWorkspace:::.isBoolGate(gh,node)){
                params <- parameters(gate)
                
                #update according to the map
                params <- sapply(params, function(param){
                      param <- gsub("<|>", "", param) #remove prefix
                      ind <- match(param, map[, "old"])
                      ifelse(is.na(ind), param, map[ind, "new"])
                    })
                names(params) <- params                          
                parameters(gate) <- params
                
              }
              negated <- flowWorkspace:::isNegated(gh, node)
                
              add(clone[[sn]], gate, name = popName, parent = thisParent, negated = negated)      
        }  
        
   }
  
  recompute(clone)
  clone
}

#' post process gs to fix channel name discrepancy caused by letter case:
#' 1.between gates defined in xml and fcs
#' 2.gates with xml
#' 
#' The assumption is the channel names is consistent across samples in FCS files.
#' So basically it loop through all gates and match up the gate parameters to 
#' the flow data and update it when needed.
fixChannelNames <- function(gs){
  
  coln <- colnames(getData(gs))
  for(sn in sampleNames(gs))
  {
    gh <- gs[[sn]]
    nodes <- getNodes(gh)
    for(node in nodes[-1])
    {
      
      g <- getGate(gh, node)
      
      if(class(g) != "booleanFilter")
      {
        isUpdate <- FALSE
        
        param <- parameters(g)
        misMatch <- which(is.na(match(param, coln)))
        #update each mismatched channel
        for(i in misMatch)
        {
          oldN <- param[i]
          matchInd <- match(tolower(oldN), tolower(coln))
          if(is.na(matchInd))
            stop(oldN, "not found in flow Data")
          else{
            isUpdate <- TRUE
            newN <- coln[matchInd]
            param[i] <- newN
            
          }
        }
        if(isUpdate){
          #update the gate
          oldParam <- parameters(g)
          
          message(paste(oldParam, col = ",")," --> ", paste(param, col = ","), "for ", node)
          parameters(g) <- param
          setGate(gh, node, g)  
        }
        
      }
      
    }
    message("done")
  }
}

#' a uility function to match fcs files based on the channels used by one example FCS file 
#' 
#' It uses \link{readFCSPar} to read parameters from FCS header to select target files, 
#' thus be used as a prefilter before \code{read.flowSet} or \link{read.ncdfFlowSet} call. 
#' 
#' @param x \code{character} vector giving the list of fcs files to match 
#' @param pattern \code{character} the example FCS file that contains the channels of interest
#' @return a \code{character} vector of fcs files that has the identical channels with \code{subset}
#' @export 
#' @examples 
#' \dontrun{
#' grep.FCS(pattern = bcells[2],  x = c(bcells,tcells))
#'  #return TRUE  TRUE FALSE FALSE
#' }
grep.FCS <- function(pattern, x){
#      browser()
      targetChnls <- readFCSPar(pattern)
      unname(sapply(x, function(thisFile){
                          thisChnls <- readFCSPar(thisFile)
                          setequal(thisChnls, targetChnls)
                          
                        })
                )
      
      
    }
#' fast way of getting channel names from fcs file by only reading header
#' 
#' This is a convenient wrapper around \link{read.FCSheader} and \code{flowCore:::readFCSgetPar}.
#' 
#' @param fileName \code{character}  fcs file name(path)
#' @return a \code{character} vector channels/parameters used in this FCS
#' @export
#' @examples 
#' 
#' readFCSPar(system.file("extdata/0877408774.B08", package = "flowCore"))
readFCSPar <- function(fileName){
  txt <- flowCore:::read.FCSheader(fileName)[[1]]
  nChannels <- as.integer(txt[["$PAR"]])
  channelNames <- unlist(lapply(1:nChannels,function(i)flowCore:::readFCSgetPar(txt,paste("$P",i,"N",sep="")))) 
  unname(channelNames)
  
}

#' a wrapper for \code{save_gslist}
#' @return a copy of original gslist with modified cdf path when cdf == "move"
#' otherwise, it behaves the same as \code{save_gslist}
#' @export 
save_gslist_labkey <- function(gslist, path, cdf, ...){
  
  save_gslist(gslist, path, cdf = cdf, ...)
  
  if(cdf == "move"){
      
      newListOfGS <- lapply(gslist, function(thisGS){
      cdfName <- basename(flowData(thisGS)@file)
      newFullPath <- file.path(path, thisGS@guid, cdfName)
      flowData(thisGS)@file <-  newFullPath
      thisGS
      }, level = 1)
      
      GatingSetList(newListOfGS)
      
      }
}

#' plot by prarent index
#' 
#' This API is mainly used for labkey module. It takes a parent index instead of the actual gate index.
#' When there is no gate associated with the x,y channel specified by user, it simply plots the \code{xyplot} 
#' or \code{densityplot} without the gate. 
#' 
#' @param x \code{character} x channel
#' @param y \code{character} y channel, if \code{NULL},then try to do \code{densityplot}
#' @export 
#' @importFrom BiocGenerics colnames

plotGate_labkey <- function(G,parentID,x,y,smooth=FALSE,cond=NULL,xlab=NULL,ylab=NULL, overlay = NULL, ...){
  #get all childrens
  cids<-getChildren(G[[1]],parentID)
  if(length(cids)>0)
  {
    #try to match to projections
#		browser()
    isMatched<-lapply(cids,function(cid){
          g<-getGate(G[[1]],cid)
          if(class(g)!="booleanFilter") 
          {
            prj<-parameters(g)
            if(length(prj)==1)#1d gate
            {
              return (prj%in%c(x,y))
              
            }else
            {
              #2d gate but y is absent
              if(is.null(y))
                return (FALSE)
              #try to match x,y to 2d gate
              revPrj<-rev(prj)
              if((x==prj[1]&&y==prj[2])||(x==revPrj[1]&&y==revPrj[2]))
                return (TRUE)
              else
                return (FALSE)	
            }
          }else
            return (FALSE)
        })
    
    ind<-which(unlist(isMatched))
    if(length(ind)>0)
      isPlotGate<-TRUE
    else
      isPlotGate<-FALSE
  }else
    isPlotGate<-FALSE
#  browser()
  formula1 <- flowWorkspace:::mkformula(c(y,x),isChar=TRUE)
#  formula1<-paste("`",y,"`~`",x,"`",sep="")
  if(!is.null(cond))
    formula1<-paste(formula1,cond,sep="|")
  formula1 <- as.formula(formula1)
#	browser()
  type <- ifelse(is.null(y), "densityplot","xyplot")
  if(isPlotGate)
    plotGate(G,cids[ind],formula=formula1,smooth=smooth,xlab=xlab,ylab=ylab, type = type, overlay = overlay, ...)
  else
  {
    fs<-getData(G,parentID)
    axisObject <- flowWorkspace:::.formatAxis(x=G[[1]],data=fs[[1]],xParam=x,yParam=y,...)
    if(is.null(xlab)){
      xlab <- axisObject$xlab
    }
    if(is.null(ylab)){
      ylab <- axisObject$ylab
    }
    if(type == "xyplot"){
      overlay <- flowWorkspace:::.getOverlay(G, overlay, params = c(x, y))
      xyplot(formula1
          ,fs
          ,smooth=smooth
          ,xlab=xlab
          ,ylab=ylab
          ,scales=axisObject$scales
          ,overlay = overlay
          ,...
              )  
    }else{
      densityplot(formula1
          ,fs
          ,xlab=xlab
          ,scales=axisObject$scales
          ,...)
    }
    
  }
  
}

#' impute the gate (flagged as failure by external algorithm) with refGate from nearest neighbour sample
#'
#' @param gs a \code{GatingSet}
#' @param node a \code{character} or \code{numeric} specifing node index
#' @param failed a \code{character} or \code{numeric} specifing samples that fails the gating QA
#' @param ... other arguments passed to \code{.nearestSample}

.nearestSamples <- function(gs, node, failed, ...){
  #get samples that do not fail the QA check
#  browser()
  samples <- sampleNames(gs)
  failedInd <- match(failed,samples)
  samples <- samples[-c(failedInd)]
  sapply(failed,function(thisTarget){
        message("Finding reference sample for: ",thisTarget)
        .nearestSample(gs, node = node, target = thisTarget, source = samples, ...) 
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


##################################################
#merge GatingSets into groups based on their gating schemes
#Be careful that the splitted resluts still points to the original data set!!
#drop the unused channels if needed before merging them 
##########################################################
.groupByTree <- function(x){
  message("Grouping by Gating tree...")
  node_seq <-unlist(lapply(x,function(this_gs){
            this_gh <- this_gs[[1]]
            this_nodes <- getNodes(this_gh, showHidden = TRUE)
            paste(this_nodes,collapse = "")
            
          }))
  split(x,node_seq)
}
#try to determine the redundant terminal nodes that can be removed
#in order to make trees mergable
.checkRedundantNodes <- function(gs_groups){
  nodeSet <- lapply(gs_groups,function(this_group){
              getNodes(this_group[[1]][[1]], showHidden = TRUE)
            })
  commonNodes <- Reduce(intersect, nodeSet)
  toRemove <- mapply(nodeSet,gs_groups,FUN=function(thisNodeSet,this_group){
                  nodesToRm <- setdiff(thisNodeSet,commonNodes)
                  #check if those nodes are terminal
                  isTerminal <- sapply(nodesToRm,function(thisNode){
                            length(getChildren(this_group[[1]][[1]],thisNode))==0
                          })
#                  browser()
                  if(!all(isTerminal)){
                    stop("Can't drop the non-terminal nodes",nodesToRm[!isTerminal])
                  }    
                  nodesToRm
                })
    toRemove       
}

#drop the terminal nodes
.dropRedundantNodes <- function(gs_groups,toRemove){
  mapply(toRemove,gs_groups,FUN=function(thisNodeSet,this_group){
        if(length(thisNodeSet)>0){
          lapply(thisNodeSet,function(thisNode){
#                browser()
                message("Removing ", thisNode)
                lapply(this_group,function(this_gs){
                      Rm(thisNode,this_gs)
                    })
              })
        }
      })
  
}

##merge gs 
#' @param force \code{logical} if TRUE, drop any channels if neccessry, 
#'                          otherwise, be conservative by only dropping unused channels
.mergeGS <- function(this_gslist, force = FALSE){
  
        
        
        
        if(length(this_gslist) > 1){
          #find the common colnames
          col_list <- lapply(this_gslist,function(this_gs)colnames(flowData(this_gs)))
          global_colnames <- Reduce(intersect, col_list)
          
          if(is.null(global_colnames))
            stop("Can't merge!no common channels.")
          
          this_gslist <- lapply(this_gslist,function(this_gs){
#                    browser()
                this_fs <- getData(this_gs)
                
                if(force){
                  toDrop <- setdiff(colnames(this_fs), global_colnames)
                  if(length(toDrop) >0)
                    message("drop ", toDrop)
                  flowData(this_gs) <- this_fs[,global_colnames]
                }else
                {
                  #drop the unused marker from fs                    
                  this_fs_colnames <- colnames(this_fs)
                  this_fr <- this_fs[[1]]
                  this_pd <- pData(parameters(this_fr))
                  within_common_chnnl <- this_fs_colnames%in%global_colnames
                  non_na_channel <- unname(!is.na(this_pd[,"desc"]))
                  to_include <- grepl(pattern="[FS]SC|[Tt]ime",this_pd[,"name"])
                  to_include <- to_include |  non_na_channel | within_common_chnnl
                  
                  if(length(which(to_include)) != nrow(this_pd)){
                    #drop channels from colnames of flowFrame
                    message("drop empty channel:",this_pd[!to_include,1])
                    fr_colnames <- colnames(this_fr)
                    fr_colnames <- fr_colnames[to_include]
                    #update the colnames of flowSet accordingly                                   
                    this_fs_colnames <- this_fs_colnames[match(fr_colnames,this_fs_colnames)]
                  }
                  
                  
                  if(!setequal(global_colnames,this_fs_colnames))
                    stop("merge failed!These channels are not common across data sets:\n"
						  , paste0(setdiff(this_fs_colnames,global_colnames))	 
                          , "\n If they are non-NA channels, try force = TRUE to drop them."
                          )
                }
                
                #reorder colnames of other gs by global_colnames
                flowData(this_gs) <- this_fs[,global_colnames] 
                
                
                this_gs
              })
        }
        
        GatingSetList(this_gslist)
      
}

#' this is the wrapper for labkey where only one gslist to be returned
#merge_gs_labkey <- function(x,...){
#  gs_groups <- .groupByTree(x, drop = TRUE)
# 
#  if(length(gs_groups)>1){
#      stop("Can't merge because multiple gating trees are present!")
#  }else{
#    res <- .mergeGS(gs_groups)
#    res[[1]]
#  }
#
#}

#TODO: to deprecate
#' cluster/merge GatingSets based on the gating tree structures.
#' 
#' merge GatingSets based on the gating tree structures.
#' 
#' @details Group the individual GatingSets by their gating schemes.It is done by comparing the node list returned by \code{\link{getNodes}},which assumes
#' they follow the same population naming conventions.
#' Meanwhile the unused channels are automatically dropped to make sure the flow data has identical data structure within each group.
#' In order to further merge multiple GatingSet objects into one, use \code{\link{rbind2}}.     
#' 
#' @param x A \code{list} of \code{GatingSet}s . 
#' @return A \code{\link{GatingSetList}} that contains multiple GatingSets each of which share the same gating and data structure.
#' @author Mike Jiang \email{wjiang2@@fhcrc.org}
#' @seealso \code{\link{rbind2}},\code{\link{GatingSetList}}
#' @examples \dontrun{
#' 	#load gatingsets from disk
#' 	#gs_toMerge is the path that stores multiple archived gatingsets
#' 	gs_list<-lapply(list.files("flowIncubator/output/gs_toMerge",full=T),function(this_folder){
#'       flowWorkspace:::load_gs(this_folder)
#'     })
#'     
#' 	gs_list <- merge_gs(gs_list)
#' 	gs_list
#' 	
#' }
#' @export 
merge_gs<-function(x,...){
#      browser()

      message("Grouping by Gating tree...")
      node_seq <-unlist(lapply(x,function(this_gs){
                this_gh <- this_gs[[1]]
                this_nodes <- getNodes(this_gh, showHidden = TRUE)
                paste(this_nodes,collapse = "")
                
              }))
      gs_groups <- split(x,node_seq)
        
      
      #start to merge      
      lapply(1:length(gs_groups),function(i){
#            browser()
            this_group <- gs_groups[[i]]

            #drop the unused marker from fs
            if(length(this_group) > 1){
              this_group <- lapply(this_group,function(this_gs){
#                    browser()
                    this_fs <- getData(this_gs)
                    
                    this_pd <- pData(parameters(this_fs[[1]]))
                    non_na_channel <- unname(!is.na(this_pd[,"desc"]))
                    to_include <- grepl(pattern="[FS]SC|[Tt]ime",this_pd[,"name"])
                    to_include <- to_include |  non_na_channel
                    if(length(which(to_include)) != nrow(this_pd)){
                    
                      message("drop empty channel:",this_pd[!to_include,1])
                      
                      flowData(this_gs) <- this_fs[,to_include]
                    
                    }
                    this_gs
                  })
            }
            
            GatingSetList(this_group)
          })
        
      
    }
      
