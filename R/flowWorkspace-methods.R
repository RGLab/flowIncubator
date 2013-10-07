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


#' routine to return the indices by specify boolean combination of reference nodes:
#' 
#' It adds the boolean gates and does the gating on the fly, and 
#' return the indices associated with that bool gate, and
#' remove the bool gate
#' the typical use case would be extracting any-cytokine-expressed cells
#' @param y a quoted expression.
#' @examples
#' getIndices(gs,quote(`4+/TNFa+|4+/IL2+`)) (it may be faster than R version)
#' @export 
#' @import flowWorkspace
setMethod("getIndices",signature=c("GatingSet","name"),function(obj, y, ...){
      
      bf <- eval(substitute(booleanFilter(v),list(v=y)))
      gh <- obj[[1]]
      
          suppressMessages({
              suppressWarnings(
                id <- add(obj,bf)
                )
                
              allNodes <- getNodes(gh,isPath=TRUE, showHidden = TRUE)
              this_node <- allNodes[id]
              
              
              res <-try(recompute(obj,id),silent=T)
           })
      
      
      if(class(res)=="try-error"){
        Rm(this_node,obj)
        stop(res)
      }else{
        this_ind <- lapply(obj,function(this_gh)getIndices(this_gh,this_node))
        Rm(this_node,obj)
        this_ind  
      }
      
   })
#' @export 
getIndiceMat<-function(gh,y){
  strExpr <- as.character(y)
  nodes <- strsplit(strExpr,split="\\|")[[1]]
#  browser()
  #extract logical indices for each cytokine gate
  indice_list <- sapply(nodes,function(this_node)getIndices(gh,this_node),simplify = FALSE)
  #construct the indice matrix
  do.call(cbind,indice_list)
}

#' create mapping between pops and channels
.getPopChnlMapping<-function(gh, y, pop_marker_list){
#  browser()
  #get pop names
  strExpr <- as.character(y)
  popNames <- strsplit(strExpr,split="\\|")[[1]]
  
  #parse the markers of interest from pop names
  markers_selected <- sapply(popNames,function(this_pop){
        this_pops <- strsplit(split="/",this_pop)[[1]]
        #get the terminal node
        term_pop <- this_pops[length(this_pops)]
        term_pop
      },USE.NAMES=FALSE)
  
  #match to the pdata of flow frame
  fr <- getData(gh, use.exprs = FALSE)
  this_pd <- pData(parameters(fr))
  all_markers <- this_pd[,"desc"]
  all_markers <- as.character(all_markers)
  all_markers[is.na(all_markers)] <- "NA"
  is_matched <- sapply(all_markers,function(this_marker){
        
        ##using the marker name provided by arguments by default
        is_manual_provided <- grep(this_marker,pop_marker_list)
        if(length(is_manual_provided)>0){
          if(length(is_manual_provided) > 1)
          {
            #more than one matched, then do the exact match
            is_manual_provided <- match(this_marker, pop_marker_list)
            if(length(is_manual_provided) > 1)
              stop(this_marker, " is matched with more than one populations")
          }
          
          res <- TRUE
          names(res) <- names(pop_marker_list)[is_manual_provided] 
        }else{
          this_matched <- grep(pattern = this_marker, x=markers_selected, fixed=TRUE)
          if(length(this_matched)>1){
            stop("multiple populations mached to:", this_marker)
          }else if(length(this_matched)==0){
            res <- FALSE
          }else{
            res <- TRUE
            names(res) <- popNames[this_matched]
          }  
        }        
          
        res
      }, USE.NAMES=FALSE)
#  browser()   
  pop_matched <- is_matched[is_matched]
  if(length(pop_matched)!=length(popNames)){
    stop("No markers in flow data matches ", "Populations:", paste(popNames[!popNames%in%names(pop_matched)],collapse="; "))
    
  }
    
  cbind(pop=names(is_matched[is_matched]),this_pd[is_matched,c("name","desc")])
  
  
  
}

#' Return the flowSet associated with a GatingSet by boolean expression
#' 
#' Returns a flowSet containing the events defined at by boolean expression \code{y}.
#' @param obj A \code{GatingSet} object .
#' @param y \code{name} boolean expression specifying the boolean combination of different cell populations
#' @return A \code{list} of \code{numerci matrices}
#' @author Mike Jiang \email{wjiang2@@fhcrc.org}
#' @seealso \code{\link{getIndices}} \code{\link{getProp}} \code{\link{getPopStats}}
#' @examples \dontrun{
#' 	#G is a GatingSet
#' 	geData(G,3)
#' 	res <- getData(gs[1],quote(`4+/TNFa+|4+/IL2+`))
#' 	res[[1]]
#' }
#' @export

setMethod("getData",signature=c("GatingSet","name"),function(obj, y,pop_marker_list = list(),...){
      #get ind of bool gate
      bool_inds <- getIndices(obj,y,...)
      
      
      lapply(obj,function(gh){
            #get pop vs channel mapping
            pop_chnl<- .getPopChnlMapping(gh,y,pop_marker_list)
            this_chnls <- as.character(pop_chnl[,"name"])
            this_pops <-  as.character(pop_chnl[,"pop"])
            
            #get mask mat
            this_sample <- getSample(gh)
            message(this_sample)
            
            this_ind <-  bool_inds[[this_sample]]
            
            if(sum(this_ind)==0){
              NULL
            }else{
              
              this_mat <- getIndiceMat(gh,y)[this_ind,this_pops, drop=FALSE]
              #subset data by channels selected
              
              this_data <- getData(gh)
              this_subset <- exprs(this_data)[this_ind,this_chnls, drop=FALSE] 
              #masking the data
              this_subset <- this_subset *  this_mat
              colnames(this_subset) <- pop_chnl[,"desc"]
              this_subset
            }
          })     
    })
setMethod("getData",signature=c("GatingSetList","name"),function(obj, y, pop_marker_list = list(), ...){
      
#      browser()
      sapply(sampleNames(obj),function(this_sample){
            message(this_sample)
            gh <- obj[[this_sample]]
            
            pop_chnl<- .getPopChnlMapping(gh,y,pop_marker_list)
            this_pops <-  as.character(pop_chnl[,"pop"])
            this_chnls <- as.character(pop_chnl[,"name"])
            
            
            #get mask mat
#      browser()
            
            this_mat <- getIndiceMat(gh,y)[,this_pops, drop=FALSE]
            #get indices of bool gates 
            this_ind <- this_mat[,1]
            for(i in 2:ncol(this_mat)){
              
              this_ind <- this_ind |this_mat[,i]
              
            }
            if(sum(this_ind)==0){
              NULL
            }else{
              this_mat <- this_mat[this_ind,,drop = FALSE]
              #subset data by channels selected
              
              this_data <- getData(gh)
              this_subset <- exprs(this_data)[this_ind,this_chnls, drop=FALSE] 
              #masking the data
              this_subset <- this_subset *  this_mat
              colnames(this_subset) <- pop_chnl[,"desc"]
              this_subset  
            }
            
          },simplify = FALSE)  
      
    })

setMethod("getIndices",signature=c("GatingSetList","name"),function(obj, y, ...){
      
    })      
##################################################
#merge GatingSets into groups based on their gating schemes
#Be careful that the splitted resluts still points to the original data set!!
#drop the unused channels if needed before merging them 
##########################################################
.groupByTree <- function(x){
  message("Grouping by Gating tree...")
  node_seq <-unlist(lapply(x,function(this_gs){
            this_gh <- this_gs[[1]]
            this_nodes <- getNodes(this_gh,isPath=T, showHidden = TRUE)
            paste(this_nodes,collapse = "")
            
          }))
  split(x,node_seq)
}
#try to determine the redundant terminal nodes that can be removed
#in order to make trees mergable
.checkRedundantNodes <- function(gs_groups){
  nodeSet <- lapply(gs_groups,function(this_group){
              getNodes(this_group[[1]][[1]],isPath=T, showHidden = TRUE)
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
                    stop("Can't merge!colnames of flow data are different. Use force = TRUE to drop any channels in order to proceed the merging")
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
                this_nodes <- getNodes(this_gh,isPath=T, showHidden = TRUE)
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
      
