#' this function retrieves the gates from GatingSet 
#' and writes a customed GatingML-2.0 file  
#' which can then be imported into cytobank  
#' operation:
#' 1. Read in gate geometry, compensation and transformation from gatingSet
#' 2. Rescale gate boundaries with flowJoTrans() so gates show up in flowJo
#' 3. Save gates and hierarchy structure to R environment
#' 4. Write environment out to gatingML using write.GatingML()
#' @import flowUtils write.gatingML
#' @import XML saveXML xmlTreeParse
#' @examples 
#' \dontrun{
#' require(flowWorkspace)
#' require(XML)
#' library(base64enc)
#' library(jsonlite)
#' localPath <- "~/rglab/workspace/openCyto"
#' gs <- load_gs(file.path(localPath,"misc/testSuite/gs-tcell_asinhtGm2"))
#' outFile <- tempfile(fileext = ".xml")
#' GatingSet2GatingML(gs, outFile)
#' }
GatingSet2GatingML <- function(gs, outFile){
  flowEnv <- GatingSet2Environment(gs) 
  tmp <- tempfile(fileext = ".xml")#ensure correct file extension for xmlTreeParse to work
  flowUtils::write.gatingML(flowEnv, tmp)
  tree <- xmlTreeParse(tmp, trim = FALSE)
  root <- xmlRoot(tree)
  # browser()
  guid_mapping <- new.env(parent = emptyenv())
  root <- addCustomInfo(root, gs, guid_mapping)
  #add pop (GateSet/BooleanAndGate)
  root <- addGateSets(root, gs, guid_mapping)  
  saveXML(root, file = outFile)
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

#' @importFrom base64enc base64encode base64decode
base64encode_cytobank <- function(x){
  x <- base64encode(charToRaw(x))
  x <- gsub("=", ".", x)
  x <- gsub("\\+", "_", x)
  x <- gsub("/", "-", x) 
  x
}
base64decode_cytobank <- function(x){
  x <- gsub("\\.", "=", x)
  x <- gsub("_", "\\+", x)
  x <- gsub("-", "/", x)
  base64decode(x)
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
      
      # trans.obj <- translist[[which(ind)]]
#       if(inverse)
#         trans.fun <- trans.obj[["inverse"]] 
#       else
#         trans.fun <- trans.obj[["transform"]] 
      #rescale
      # gate <- transform(gate, trans.fun, param)
      
      #attach trans reference
      gate@parameters[[i]] <- flowEnv[[transNames[ind]]]
    }else if(nMatched > 1)
      stop("multiple trans matched to :", param)
  }
  
  gate
  
}
#' @import XML xmlTree
addGateSets <- function(root, gs, ...)
{
  
  nodePaths <- getNodes(gs, showHidden = TRUE)[-1]
  # browser()
  newNodes <- lapply(seq_along(nodePaths), function(gate_id){
                      nodePath <- nodePaths[gate_id]
                      curNode <- nodePath
                      pop_name <- basename(nodePath)
                      gate_id_path <- gate_id
                      # browser()
                      repeat{
                        curNode <- getParent(gs, curNode)
                        if(curNode == "root")
                          break
                        else{
                          cur_parent_id <- match(curNode, nodePaths)
                          gate_id_path <- c(cur_parent_id, gate_id_path)
                        }
                          
                      }
                      GateSetNode(gate_id, pop_name, gate_id_path, nodePaths, ...)
                    })
  
  addChildren(root, kids = newNodes)
}

#' @importFrom jsonlite toJSON
GateSetNode <- function(gate_id, pop_name, gate_id_path, nodePaths, guid_mapping){

  attrs = c("gating:id" = paste("GateSet", gate_id, sep = "_"))
  
  definition <- toJSON(list(gates = gate_id_path, negGates = vector()))
  
  #duplicate the refs if it is the root
  ref_gate_id_path <- gate_id_path
  if(length(ref_gate_id_path) == 1)
    ref_gate_id_path <- c(ref_gate_id_path, ref_gate_id_path)
  xmlNode("gating:BooleanGate", attrs = attrs
          , xmlNode("data-type:custom_info"
                    , xmlNode("cytobank"
                              , xmlNode("name", pop_name)
                              , xmlNode("gate_set_id", gate_id)
                              , xmlNode("definition", I(definition))#set AsIs to avoid xml escaping
                              )
                    )
          
         ,  xmlNode("gating:and"
                  #create two dummy reference
                  , .children = lapply(ref_gate_id_path, function(gate_id){
                    
                    guid <- guid_mapping[[as.character(gate_id)]]
                    attrs = c("gating:ref" = guid)
                    xmlNode("gating:gateReference", attrs = attrs)  
                  })
                )
        )
}

#' add customInfo nodes to each gate node and add BooleanAndGates
#' @import XML xmlAttrs getNodeSet
addCustomInfo <- function(root, gs, guid_mapping){
  nodePaths <- getNodes(gs, showHidden = TRUE)[-1]
  fcs_names <- pData(gs)[["name"]]

  
  for(id in 1:length(root)){
    
    curNode <- root[[id]]
    guid <- as.vector(xmlAttrs(curNode, "gating:id"))
    if(!is.null(guid)&&grepl("gate_", guid)){
        # browser()
        fields <- strsplit(guid, "_")[[1]]
        gate_id <- as.integer(fields[[2]])
        fcs_id <- as.integer(fields[[3]])
        
        nodePath <- nodePaths[gate_id]
        pop_name<- basename(nodePath)
        fcs_name <- ifelse(fcs_id == 1, "", fcs_names[fcs_id])
        gate <- getGate(gs[[1]], nodePath)
        gate_type <- class(gate)
        if(gate_type == "rectangleGate"){
          if(length(parameters(gate)) == 1)
            gate_type <- "RangeGate"
          else
            gate_type <- "RectangleGate"
        }else if(gate_type == "polygonGate")
          gate_type <- "PolygonGate"
        else if(gate_type == "ellipsoidGate")
          gate_type <- "EllipseGate"
        else
          stop("unsupported gate: ", gate_type)
        # browser()
        #insert custom info
        customNode <- customInfoNodeForGate(id, gate_id, pop_name, fcs_name, gate_type)
        newNode <- addChildren(curNode, kids = list(customNode), at = 0)        
        #modify id
        guid.new <- paste("Gate", id, base64encode_cytobank(pop_name), sep = "_")
        xmlAttrs(newNode, suppressNamespaceWarning = TRUE, append = FALSE) <- c("gating:id" = guid.new)
        #update the tree
        root[[id]] <- newNode  
        
        #record the mapping between gate_id and guid.new for the refs of GateSets
        if(fcs_id == 1)
          guid_mapping[[as.character(gate_id)]] <- guid.new
    }
  }
  root
  
}

#' @import XML newXMLNode
customInfoNodeForGate <- function(id, gate_id, pop_name, fcs_name, type)
{
  myscale <- list(flag = 1
                  , argument = "1"
                  , min = 0 
                  , max = 260000.0
                  , bins = 256
                  , size = 256
  )
  definition <- list(scale = list())
  if(type == "RangeGate")
    definition[["scale"]] <- myscale
  else{
    definition[["scale"]][["x"]] <- myscale
    definition[["scale"]][["y"]] <- myscale
  }
    
  definition <- toJSON(definition, auto_unbox = TRUE)
 #avoid using newXMLNode since it is not segfault-free.
  xmlNode("data-type:custom_info"
      , xmlNode("cytobank"
          , xmlNode("name", pop_name)
          , xmlNode("id", id)
          , xmlNode("gate_id", gate_id)
          , xmlNode("type", type)
          , xmlNode("fcs_file_filename", fcs_name)
          , xmlNode("definition", I(definition))
          )
      )
}

