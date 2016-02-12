#' this function retrieves the gates from GatingSet 
#' and writes a customed GatingML-2.0 file  
#' which can then be imported into cytobank  
#' operation:
#' 1. Read in gate geometry, compensation and transformation from gatingSet
#' 2. Rescale gate boundaries with flowJoTrans() so gates show up in flowJo
#' 3. Save gates and hierarchy structure to R environment
#' 4. Write environment out to gatingML using write.GatingML()
#' @examples 
#' \dontrun{
#' require(flowWorkspace)
#' require(XML)
#' localPath <- "~/rglab/workspace/openCyto"
#' gs <- load_gs(file.path(localPath,"misc/testSuite/gs-tcell_asinhtGm2"))
#' outFile <- tempfile(fileext = ".xml")
#' GatingSet2GatingML(gs, outFile)
#' }
GatingSet2GatingML <- function(gs, outFile){
  flowEnv <- GatingSet2Environment(gs) 
  tmp <- tempfile(fileext = ".xml")#ensure correct file extension for xmlTreeParse to work
  flowUtils::write.gatingML(flowEnv, tmp)
  tree <- xmlTreeParse(tmp, useInternalNodes = TRUE)
  addCustomInfo(tree, gs)
  #add pop (GateSet/BooleanAndGate)
  addGateSet(tree, gs)  
  saveXML(tree, file = outFile)
}


#' extract the trans, comps and gates into GML2 compatible format
#' and save to environment
GatingSet2Environment <- function(gs) {
  flowEnv <- new.env(parent = emptyenv())
  
  #parse comp and channel names
  comp <- gs@compensation
  if(is.null(comp)){
    #(assume it is identical across samples)
    comp.obj <- flowWorkspace:::.cpp_getCompensation(gs@pointer, sampleNames(gs)[1])  
    if(is.null(comp.obj)){
      #no compensation and channel names in transformation are not prefixed
      chnls <- colnames(getData(gs))
    }else{
      #parsed from flowJo and channel names are usually prefixed
      #thus get the raw channel names from here
      chnls <- comp.obj[["parameters"]]
      comp <- compensation(matrix(comp.obj$spillOver
                                  ,nrow=length(chnls)
                                  ,ncol=length(chnls)
                                  ,byrow=TRUE
                                  ,dimnames=list(chnls,chnls)
                                  )
                           )
      
    }
  }else{
    #compensation was added in R
    #channel names are not prefixed
    chnls <- as.vector(parameters(comp))
  }
  compId <- identifier(comp)
  compId <- paste("Spill", compId, sep = "_")
  identifier(comp) <- compId
  #add comp
  if(!is.null(comp))
    flowEnv[[compId]] <- comp
  
  #add trans (assume it is identical across samples)
  trans <- getTransformations(gs[[1]], only.function = FALSE)
  for(transName in names(trans)){
    trans.obj <- trans[[transName]]
    type <- trans.obj[["name"]]
    ind <- sapply(chnls, grepl, x = transName, USE.NAMES = FALSE)
    nMatched <- sum(ind)
    if(nMatched > 1)
      stop("More than one channels matched to transformation: ", transName)
    else if(nMatched == 1){
      chnl <- chnls[ind]
      if(type == "asinhtGml2"){
        #extract parameters
        trans.func <- trans.obj[["transform"]]
        env <- environment(trans.func)
        asinhtGml2.obj <- asinhtGml2(parameters = compensatedParameter(chnl
                                                                       , spillRefId = compId
                                                                       , searchEnv = flowEnv)
                                     , M = env[["m"]]
                                     , T = env[["t"]]
                                     , A = env[["a"]]
                                     , transformationId = paste0("Tr_Arcsinh_", transName)
                                    )
        
        flowEnv[[transName]] <- asinhtGml2.obj
      }else
        stop("unsupported trans: ", type)
      
    }
    
  }
  
  #add gates and pops(as GateSets)
  nodePaths <- getNodes(gs, showHidden = TRUE)[-1]
  fcs_names <- pData(gs)[["name"]]
  for(gate_id in seq_along(nodePaths)){
    nodePath <- nodePaths[gate_id]
    gates <- getGate(gs, nodePath)
    popName <- basename(nodePath)
    for(fcs_id in seq_along(fcs_names)){
      sn <- fcs_names[fcs_id]
      gate <- gates[[sn]]
# browser()      
      #transform to raw scale
      #and attach comp and trans reference to parameters
      gate <- processGate(gate, trans, inverse = TRUE, flowEnv)

      parent <- getParent(gs, nodePath)
      if(parent == "root")
        parent_id <- 0
      else
        parent_id <- match(parent, nodePaths)
      
      guid <- paste("gate", gate_id, fcs_id, sep = "_")
      identifier(gate) <- guid
      #add gate
      flowEnv[[guid]] <- gate
      
    }
  }
  flowEnv
}
setMethod("transform", signature = c("polygonGate"), function(`_data`, ...){
  .transform.polygonGate(`_data`, ...)
})
.transform.polygonGate <- function(gate, trans.fun, param){
  browser()
}

setMethod("transform", signature = c("rectangleGate"), function(`_data`, ...){
  .transform.rectangleGate(`_data`, ...)
})
.transform.rectangleGate <- function(gate, trans.fun, param){
  
      min <- gate@min[[param]]
      if(!is.infinite(min))
        gate@min[[param]] <- trans.fun(min)
      
      max <- gate@max[[param]]
      if(!is.infinite(max))
        gate@max[[param]] <- trans.fun(max)
  
  gate
  
}
processGate <- function(gate, translist, inverse = FALSE, flowEnv){
  
  params <- as.vector(parameters(gate))
  transNames <- names(translist)
  chnls <- transNames
  for(i in seq_along(params)){
    param <- params[i]
    ind <- sapply(chnls, function(chnl)grepl(chnl, param), USE.NAMES = FALSE)
    nMatched <- sum(ind)
    if(nMatched == 1){
      
      trans.obj <- translist[[which(ind)]]
      if(inverse)
        trans.fun <- trans.obj[["inverse"]] 
      else
        trans.fun <- trans.obj[["transform"]] 
      #rescale
      gate <- transform(gate, trans.fun, param)
      
      #attach trans reference
      transID <- 
      gate@parameters[[i]] <- flowEnv[[transNames[ind]]]
    }else if(nMatched > 1)
      stop("multiple trans matched to :", param)
  }
  
  gate
  
}
addGateSet <- function(doc, gs)
{
  tree <- xmlTree(doc = doc)
  nodePaths <- getNodes(gs, showHidden = TRUE)[-1]
  for(gate_id in seq_along(nodePaths)){
    nodePath <- nodePaths[gate_id]
    curNode <- nodePath
    pop_name <- basename(nodePath)
    gate_id_path <- as.character(gate_id)
    # browser()
    repeat{
      curNode <- getParent(gs, curNode)
      if(curNode == "root")
        break
      else{
        cur_parent_id <- match(curNode, nodePaths)
        gate_id_path <- paste(cur_parent_id, gate_id_path, sep = ",")
      }
        
    }
    
    attrs = c("gating:id" = paste("GateSet", gate_id, sep = "_"))
    # browser()
    tree$addNode("gating:BooleanGate", attrs = attrs, close = FALSE)
    tree$addNode("data-type:custom_info", close = FALSE)
    tree$addNode("cytobank", close = FALSE)
    tree$addNode("name", pop_name)
    tree$addNode("gate_set_id", gate_id)
    tree$addNode("definition", paste0("{'gates':[", gate_id_path, "],'negGates':[]}"))
    tree$closeTag()
    tree$closeTag()
    tree$addNode("gating:and", close = FALSE)
    #create two dummy reference
    for(i in 1:2){
      attrs = c("gating:ref" = paste("gate", gate_id, "1", sep = "_"))
      tree$addNode("gating:gateReference", attrs = attrs)  
    }
    tree$closeTag()
    tree$closeTag()
  }
}


#' add customInfo nodes to each gate node and add BooleanAndGates
addCustomInfo <- function(tree, gs){
  nodePaths <- getNodes(gs, showHidden = TRUE)[-1]
  fcs_names <- pData(gs)[["name"]]
    
  gateNodes <- getNodeSet(tree, path = paste0("/gating:Gating-ML/*[contains(@gating:id,'gate_')]"))
  for(gateNode in gateNodes){
    # browser()
    guid <- as.vector(xmlAttrs(gateNode, "gating:id"))
    fields <- strsplit(guid, "_")[[1]]
    gate_id <- as.integer(fields[[2]])
    fcs_id <- as.integer(fields[[3]])
    
    nodePath <- nodePaths[gate_id]
    pop_name<- basename(nodePath)
    fcs_name <- ifelse(fcs_id == 1, "", fcs_names[fcs_id])
    customInfo <- customInfoNodeForGate(gate_id, pop_name, fcs_name)
    addChildren(gateNode, kids = list(customInfo), at = 0)
  }
    
  
}
customInfoNodeForGate <- function(gate_id, pop_name, fcs_name){
  newXMLNode("data-type:custom_info"
          , newXMLNode("cytobank"
                    , newXMLNode("name", pop_name)
                    , newXMLNode("gate_id", gate_id)
                    , newXMLNode("fcs_file_filename", fcs_name)
                    , newXMLNode("definition", "{'scale':{'x':{'flag':4,'argument':'5','min':-20,'max':10000.0,'bins':256,'size':256},'y':{'flag':4,'argument':'5','min':-20.0,'max':10000.0,'bins':256,'size':256}}}")
  
                  )
        )
}

