#' It is essentially a graphNEL class and exists for the purpose of method dispatching.
#' @exportClass graphGatingML
setClass("graphGatingML", contains = "graphNEL")

#' Parser for gatingML exported by Cytobank
#' 
#' The Default parser (flowUtils::read.gatingML) does not  parse the population tree as well as 
#' the custom information from cytobank. (e.g. gate name, fcs filename).
#' 
#' @param file Gating-ML XML file
#' @param ... additional arguments passed to the handlers of 'xmlTreeParse'
#' @export
#' @importFrom flowUtils read.gatingML
#' @return a graphNEL that represents the population tree. 
#' The gate and population name are stored in nodeData of each node. 
#' Compensation and transformations are stored in graphData.
#' @examples 
#' 
#' xml <- system.file("extdata/cytotrol_tcell_cytobank.xml", package = "flowIncubator")
#' g <- read.gatingML.cytobank(xml) #parse the population tree
#' plotTree(g) #visualize it
#' nodeData(g)[[1]] # access individual gates
#' g@@graphData[["compensation"]] #access comps
#' g@@graphData[["transformations"]] #access trans
read.gatingML.cytobank <- function(file, ...){
  
  #parse all the elements:gate, GateSets, comp, trans
  flowEnv <- new.env()
  read.gatingML(file, flowEnv) 
  
  #parse gate info (id vs fcs and pop name)
  gateInfo <- parse.gateInfo(file)
  
  #construct tree from GateSets
  g <- constructTree(flowEnv, gateInfo)
  
  #attach comp and trans
  
  objNames <- ls(flowEnv)
  
  trans <-  sapply(objNames, function(i){
                        
                              obj <- flowEnv[[i]]
                              if(extends(class(obj), "transformation"))
                                obj
                            }, USE.NAMES = FALSE)
  g@graphData[["transformations"]] <- compact(trans)
  
  comps <- sapply(objNames, function(i){
    
                                        obj <- flowEnv[[i]]
                                        if(class(obj) == "compensation")
                                          obj
                                      }, USE.NAMES = FALSE)
  g@graphData[["compensation"]] <- compact(comps)
  
  as(g, "graphGatingML")
  
}

#' Parse the cytobank custom_info for each gate
#' 
#' Fcs filename and gate name stored in 'custom_info' element are beyong the scope of
#' the gatingML standard and thus not covered by the default 'read.gatingML'. 
#' 
#' @param file xml file path
#' @param ... additional arguments passed to the handlers of 'xmlTreeParse'
#' @return a data.frame that contains three columns: id (gateId), name (gate name), fcs (fcs_file_filename).
#' @importFrom XML xmlRoot xmlName xmlGetAttr xmlValue xmlElementsByTagName xmlChildren
#' @importFrom plyr ldply
parse.gateInfo <- function(file, ...)
{       
  
  root <- xmlRoot(flowUtils:::smartTreeParse(file,...))
  
  ldply(xmlChildren(root), function(node){
    nodeName <- xmlName(node)
    
    if(grepl("*Gate", nodeName)){
      # if(nodeName=="BooleanGate"){
        id <- xmlGetAttr(node, "id")
        
        name <- getCustomNodeInfo(node, "name")
        fcs_file_filename <- getCustomNodeInfo(node, "fcs_file_filename")
        gate_id <- getCustomNodeInfo(node, "gate_id")
        
        
        c(id = id, name = name, gate_id = gate_id, fcs = fcs_file_filename)
        
      # }
    }
    
  }, .id = NULL)
}

getCustomNodeInfo <- function(node, nodeName){
  custom_node <- xmlElementsByTagName(node, nodeName, recursive = TRUE)
  if(length(custom_node) >0 ){
    value <- xmlValue(custom_node[[1]])
    if(length(value)==0)
      value <- ""
  }
  else 
    value <- ""
  
  value
}

#' Given the leaf node, try to find out if a collection of nodes can be matched to a path in a graph(tree) by the bottom-up searching
#' @param g graphNEL
#' @param leaf the name of leaf(terminal) node 
#' @param nodeSet a set of node names
#' @return TRUE if path is found, FALSE if not path is matched.
#' @importFrom graph inEdges
matchPath <- function(g, leaf, nodeSet){
  
    if(leaf %in% nodeSet){
      nodeSet <- nodeSet[-match(leaf, nodeSet)] #pop it out 
      leaf <- inEdges(leaf, g)[[1]]#update it with parent node
      if(length(nodeSet) == 0){ #nodeSet emptied
        if(length(leaf) == 0)  #and path reach to the top
          return (TRUE)
        else
          return(FALSE) #path has not reached to the top
      }else
      {
        if(length(leaf) == 0) #reach the top before nodeSet is emptied
          return(FALSE)
        else
          matchPath(g, leaf, nodeSet)  #continue match the path against nodeSet  
      }
    }else
      return(FALSE)

}

#' Reconstruct the population tree from the GateSets
#' @param flowEnv the enivornment contains the elements parsed by read.gatingML function
#' @param gateInfo the data.frame contains the gate name, fcs filename parsed by parse.gateInfo function
#' @return a graphNEL represent the population tree. The gate and population name are stored as nodeData in each node.
#' @importFrom graph graphNEL nodeDataDefaults<- nodeData<- addEdge edges removeNode
#' @importFrom plyr compact dlply
constructTree <- function(flowEnv, gateInfo){
  
  #parse the references from boolean gates
  #currently we use refs (which point to the gate id defined by the standard)
  #we may want to switch to the gate_id in custom_info used by cytobank internally
  gateSets <- sapply(ls(flowEnv), function(i){
    obj <- flowEnv[[i]]
    if(class(obj) == "intersectFilter"){
      refs <- obj@filters
      refs <- sapply(refs, slot, "name")
      unique(refs) #make root node depth of 1
    }
  })
  
  gateSets <- compact(gateSets)
  popIds <- names(gateSets)
  #sort by the depths of path
  counts <- sapply(gateSets, function(i)length(i))
  counts <- sort(counts)
  
  #init the graph with all gates that have been referred
  # gateIds <- unique(unlist(gateSets))
  g <- graphNEL(nodes = popIds, edgemode = "directed")
  
  nodeDataDefaults(g, "popName") <- ""
  nodeDataDefaults(g, "gateInfo") <- list()
  
  # add edges (from nDepths 2 to n)
  for(popId in names(counts)){
    thisGateSet <- gateSets[[popId]]
    nDepth <- length(thisGateSet)
    
    popName <- subset(gateInfo, id == popId)[["name"]]
    #add pop name
    nodeData(g, popId, "popName") <- popName
    
    if(nDepth == 1){#root node
      gateID <- thisGateSet[1]

    }else{
      #traverse all the potential parent GateSets(i.e. subset)
      parents_popIds <- names(counts[counts == nDepth - 1])
      for(parentID in parents_popIds)
      {
        pGateSet <- gateSets[[parentID]]
        #check the overlaps between two sets
        if(setequal(intersect(pGateSet, thisGateSet), pGateSet))
        { #if it is indeed a subset of the current gateSet
          gateID <- setdiff(thisGateSet, pGateSet) #parse out the gateID
          #add edge
          g <- addEdge(parentID, popId, g)
          # There might be multiple parents due to the duplicated GateSets allowed in gatingML
          # it is not clear which one should be picked and
          # we pick the first one for now
          break 
        }
      }
    }
    
    #add gate
    sb <- subset(gateInfo, id == gateID)
    #try to find the tailored gate
    tg_sb <- subset(gateInfo, gate_id == sb[["gate_id"]] & fcs != "")
    
    tg <- dlply(tg_sb, "id", function(row){
                gateID <- row[["id"]]
                flowEnv[[gateID]]
              })
    names(tg) <- tg_sb[["fcs"]]
    nodeData(g, popId, "gateInfo") <- list(list(gate = flowEnv[[gateID]]
                                                , gateName = sb[["name"]]
                                                , fcs = sb[["fcs"]]
                                                , tailored_gate = tg
                                              )
                                          )
    
  }
  
  #prune the tree by removing the orphan nodes
#   egs <- edges(g)
#   leaves <- names(egs[sapply(egs, length) == 0])
#   for(node in leaves){
#     if(length(inEdges(node, g)[[1]]) == 0)
#       g <- removeNode(node, g)
#   }
  g
}

#' plot the population tree
#' @param x a graphNEL generated by constructTree function
#' @param y not used
#' @param label specifies what to be dispaled as node label. Can be either 'popName' (population name parsed from GateSets) or 'gateName'(the name of the actual gate associated with each node)
#' @export
#' @importFrom graph nodeData nodes<- nodeRenderInfo<-
#' @importFrom Rgraphviz renderGraph layoutGraph
setMethod("plot", signature = c(x = "graphGatingML", y = "missing"), definition = function(x, label = c("popName", "gateName")){
  label <- match.arg(label, c("popName", "gateName"))
  if(label == "popName")
    nodeLabel  <- sapply(nodeData(x), `[[`, "popName")
  else
    nodeLabel  <- sapply(nodeData(x), function(i)i[["gateInfo"]][["gateName"]])
  
  
  #annotate the node with tailor gate info
  nTailoredGate <- sapply(nodeData(x), function(i)length(i[["gateInfo"]][["tailored_gate"]]))
                    
  nAttrs <- list()
  
  nAttrs$label <- nodeLabel

  nAttrs$lty <- sapply(nTailoredGate
                       ,function(i)
                       {
                         ifelse(i>0,"dotted","solid")
                       })
  
  nodeRenderInfo(x) <- nAttrs
  lay <- layoutGraph(x
                     ,attrs=list(graph=list(rankdir="LR",page=c(8.5,11))
                             ,node=list(fixedsize=FALSE
                                      ,fontsize = 12
                                      ,shape="ellipse"
                           )
               )
              )
  renderGraph(lay)
  
})

#' Apply the gatingML graph to a GatingSet
#' It performs compensation, transformation and gating.
#' @param x graphGatingML
#' @param y GatingSet
#' @param ... other arguments
setMethod("gating", c("graphGatingML", "GatingSet"), function(x, y, ...){
  gating.graphGatingML(x, y, ...)
})

gating.graphGatingML <- function(x, y, ...) {
  
  fs <- flowData(gs)
  message("compensation...")
  comp <- x@graphData[["compensation"]][[1]]
  ##TODO: determine the comp based on compensation-ref
  param <- parameters(comp)
  param <- gsub("Comp_", "", param)
  colnames(comp@spillover) <- param
  fs <- compensate(fs, comp)
  
  message("transformation...")
  trans <- x@graphData[["transformations"]]
  #remove the generic one
  samples <- sampleNames(fs)
  for(tran in trans){
    
    chnl <- as.vector(parameters(tran@parameters))
    if(chnl!="any"){
      for(sn in samples){
        fr <- fs[[sn]]
        exprs(fs[[sn]])[,chnl] <- eval(tran)(fr)
      }
    }
    
  }
  flowData(gs) <- fs
  #can't use transform method 
  #since Gml version of trans object does not follow the convention of the transformation function (input)
  #thus can't construct the valid transformList to be operated on
  # chnl <- sapply(trans, function(tran)unname(parameters(tran@parameters)), USE.NAMES = F)
  # trans <- transformList(chnl,trans)
  # fs <- transform(fs, trans)
  
  
  message("gating...")
  gs
}
