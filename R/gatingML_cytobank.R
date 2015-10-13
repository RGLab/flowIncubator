
# parse.sampleName <- function(file, ...)
# {       
#   require(XML)
#   require(plyr)
#   root <- xmlRoot(flowUtils:::smartTreeParse(file,...))
#   
#   ldply(xmlChildren(root), function(node){
#     nodeName <- xmlName(node)
#     
#     if(grepl("*Gate", nodeName)){
#       if(nodeName!="BooleanGate"){
#         gateId <- xmlGetAttr(node, "id")
#         
#         fcs_name <- xmlValue(xmlElementsByTagName(node,"fcs_file_filename", recursive = TRUE)[[1]])
#         
#         if(length(fcs_name)>0)
#           c(id = gateId, fcs = fcs_name)
#           
#       }
#     }
#     
#   }, .id = NULL)
# }
#-------#default read.gatingML doesn't parse fcs and name info -----------------------
parse.gateInfo <- function(file, ...)
{       
  require(XML)
  require(plyr)
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


matchPath <- function(g, leaf, queue){
  
    if(leaf %in% queue){
      queue <- queue[-match(leaf, queue)] #pop it out 
      leaf <- inEdges(leaf, g)[[1]]#update it with parent node
      if(length(queue) == 0){ #queue emptied
        if(length(leaf) == 0)  #and path reach to the top
          return (TRUE)
        else
          return(FALSE) #path has not reached to the top
      }else
      {
        if(length(leaf) == 0) #reach the top before queue is emptied
          return(FALSE)
        else
          matchPath(g, leaf, queue)  #continue match the path against queue  
      }
    }else
      return(FALSE)

}

#discover trees from parsed gatingML gates
constructTree <- function(flowEnv, gateInfo){
  objs <- as.list(flowEnv) #convert to list for easy operation
  #get boolean gates
  pops <- lapply(objs, function(obj){
    if(class(obj) == "intersectFilter"){
      refs <- obj@filters
      refs <- sapply(refs, slot, "name")
      refs
    }
  })
  pops <- pops[!sapply(pops,is.null)]
  
  #sort by nCombinations
  counts <- sapply(pops, length)
  counts <- sort(counts)
  
  #init the graph with all refered gates
  allNodes <- unique(unlist(pops))
  g <- graphNEL(nodes = allNodes, edgemode = "directed")
  
  nodeDataDefaults(g, "popName") <- ""
  nodeDataDefaults(g, "gateInfo") <- list()
  
  # add edges
  for(popId in names(counts)){
    thisRef <- pops[[popId]]
    thisCount <- length(thisRef)
    
    popName <- subset(gateInfo, id == popId)[["name"]]

    if(thisCount == 2){
      #assuming thisRef[1] has already been (or will be) taken care of
      #which means we assume root node always was defined by single Gate
      nodeId <- thisRef[2]
      #add gate info
      sb <- subset(gateInfo, id == nodeId)
      nodeData(g, nodeId, "gateInfo") <- list(list(gate = flowEnv[[nodeId]], gateName = sb[["name"]], fcs = sb[["fcs"]]))
      
      #add pop Info
      nodeData(g, nodeId, "popName") <- popName
      
      if(thisRef[1] != nodeId)#root node defined by single gate
      {
        g <- addEdge(thisRef[1], nodeId, g)
      }
      
      
    }else{
      #traverse all the potential parent GateSets(i.e. subset)
      parents <- names(counts[counts == thisCount - 1])
      for(parent in parents)
      {
        pRef <- refs[[parent]]
        if(length(unique(pRef))>1)
          if(setequal(intersect(pRef, thisRef), pRef))
          { #if it is indeed a subset of the current gateSet
            nodeId <- setdiff(thisRef, pRef) #parse out the leave node first
            #then look up the existing graph edges to find its direct parent node
            for(leaf in edges(g))
            {
              
              if(length(leaf)>0){
                if(matchPath(g, leaf, pRef)){
                  #add gate info 
                  sb <- subset(gateInfo, id == nodeId)
                  nodeData(g, nodeId, "gateInfo") <- list(list(gate = flowEnv[[nodeId]], gateName = sb[["name"]], fcs = sb[["fcs"]]))
                  #add pop Info
                  nodeData(g, nodeId, "popName") <- popName
                  #add edge
                  g <- addEdge(leaf, nodeId, g)
                }
                  
                
              }
            }
            
          }
      }
      
    }
    
  }
  
  #prune the tree by removing the orphan nodes
  egs <- edges(g)
  leaves <- names(egs[sapply(egs, length) == 0])
  for(node in leaves){
    if(length(inEdges(node, g)[[1]]) == 0)
      g <- removeNode(node, g)
  }
  g
}

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
