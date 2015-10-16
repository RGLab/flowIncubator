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
  
  g
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
        gateId <- xmlGetAttr(node, "id")
        
        pop <- xmlValue(xmlElementsByTagName(node,"name", recursive = TRUE)[[1]])
        
        fcsNodes <- xmlElementsByTagName(node,"fcs_file_filename", recursive = TRUE)
        if(length(fcsNodes) >0 ){
          fcs <- xmlValue(fcsNodes[[1]])
          if(length(fcs)==0)
            fcs <- ""
        }
        else 
          fcs <- ""
    
      c(id = gateId, name = pop, fcs = fcs)
        
      # }
    }
    
  }, .id = NULL)
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
#' @importFrom plyr compact
constructTree <- function(flowEnv, gateInfo){
  
  #get boolean gates
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
    nodeData(g, popId, "gateInfo") <- list(list(gate = flowEnv[[gateID]]
                                                , gateName = sb[["name"]]
                                                , fcs = sb[["fcs"]])
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
#' @param g a graphNEL generated by constructTree function
#' @param label specifies what to be dispaled as node label. Can be either 'popName' (population name parsed from GateSets) or 'gateName'(the name of the actual gate associated with each node)
#' @export
#' @importFrom graph nodeData nodes<-
plotTree <- function(g, label = c("popName", "gateName")){
  label <- match.arg(label, c("popName", "gateName"))
  if(label == "popName")
    nodeLabel  <- sapply(nodeData(g), `[[`, "popName")
  else
    nodeLabel  <- sapply(nodeData(g), function(x)x[["gateInfo"]][["gateName"]])
  nodes(g) <- nodeLabel
  plot(g, attrs=list(node=list(shape = "ellipse"
                               , fontsize = 12
                               , fixedsize = F
                               # , label = "nodeLabel" #somehow label attr is not working
                               )
                     , graph=list(rankdir="LR")))
  
}

