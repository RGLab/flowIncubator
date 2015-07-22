gsGatingML <- function (gs, outFile) {

# this function retrieves flowWorkspace gatingSet 
# and writes a valid GatingML-2.0 file using flowUtils package 
# which can then be imported into flowJo  

# operation:
# 1. Read in gate geometry, compensation and transformation from gatingSet
# 2. Rescale gate boundaries with flowJoTrans() so gates show up in flowJo
# 3. Save gates and hierarchy structure to R environment
# 4. Write environment out to gatingML using write.GatingML()

# Phu T. Van, FHCRC, w/ lots of help from G. Finak
# and input from J. Spidlen, J. Frelinger and J. Almarode

# load relevant libraries
require(flowWorkspace)
require(flowUtils)

# check that the gatingSet is valid
if(!is(gs, "GatingSet")) { 
  stop ("Invalid gatingSet supplied !!! STOPPING...")
} else {
  message("gatingSet appears to be valid, parsing...\n")
  
}
  
# get gates from our gatingSet, exclude `root`
gates <- lapply(getNodes(gs)[-1], function(x)getGate(gs[[1]], x))
gateNames <- getNodes(gs)[-1]
names(gates) <- gateNames

# get each gate's parent, so we can reconstruct the hierarchy later
parents <- lapply(names(gates), function(x)getParent(gs, x))
names(parents) <- gateNames

# CLEAN GATE NAMES 
# (remove <> from marker names)
for (i in 1:length(gates)) { 
  s <- gates[[i]]
  parameters(s) <- gsub("<|>", "", parameters(s))
  names(parameters(s)) <- gsub("<|>", "", names(parameters(s)))
  
  if (is(s, "polygonGate")){
    colnames(s@boundaries) <- gsub("<|>", "", colnames(s@boundaries))
  
  } else if (is(s, "rectangleGate")){
    names(s@min) <- gsub("<|>", "", names(s@min))
    names(s@max) <- gsub("<|>", "", names(s@max))
    
  } else if (is(s, "ellipsoidGate")){
    names(s@mean) <- gsub("<|>", "", names(s@mean))
    colnames(s@cov) <- gsub("<|>", "", colnames(s@cov))
    rownames(s@cov) <- gsub("<|>", "", rownames(s@cov))
  }
  gates[[i]] <- s
}

# GET MARKER INFORMATION
# use.exprs = F gives us markers w/o having to read CDF 
# keep track of which channels need to be transformed later
message("getting marker information...\n ")
markers <- pData(parameters(getData(gs[[1]],  use.exprs = F)))
rownames(markers) <- markers$name
rownames(markers) <- gsub("<|>", "", rownames(markers))
markers <- markers[,!(names(markers) %in% "name")]
scaleme <- rownames(markers[!is.na(markers["desc"]),])


# GET COMPENSATION MATRICES FROM GATINGSET
spillMat <- getCompensationMatrices(gs[[1]])

if (!is.null(spillMat)){
    spillMat <- as.matrix(spillMat@spillover)
} else { # no compensation matrix, eg cyTOF, create identity matrix
    spillMat <- matrix(data=1.0, nrow = length(scaleme), 
                                 ncol = length(scaleme))
    colnames(spillMat) <- sort(scaleme)
    rownames(spillMat) <- sort(scaleme)
}

compMat <- compensation(spillover=spillMat, compensationId='compMat')

# get transformation functions from gs, clean channel/marker names
# ie. replace spaces with "_", prepend "tr" 
transforms <- getTransformations(gs[[1]])

names(transforms) <- unlist(regmatches(names(transforms), 
                            gregexpr("(<.*?>)", names(transforms)))
                            )
names(transforms) <- gsub("<|>", "", names(transforms))

# make the environment for flowUtils, assign the transforms 
flowEnv <- new.env()

# remake transforms using gatingML's logicleGml2 and assign them to the env
# logicle parameters estimated from flowJoTrans by G.Finak
# W = 0.3176969  linear decades, ~0.5 for flowJo
# M = 4.4941027  log decades,   = 4.5 for flowJo
# A = 0.2653686  negative decades,   0
# T = 262144     maximum data value

# NOTE: as of flowUtils 1.30 if you supply multiple different transformations with the  
# same parameterization, flowUtils will only define ONE transform (the first one), then
# refer to it using gating:transformation-ref="YourFirstTransform", ignoring the 
# remaining transformations. This is Intended Behavior for GatingML 2.0

for (i in 1:length(transforms)){
  t <- logicletGml2(parameters = names(transforms)[[i]],
                        W = 0.5,  
                        M = 4.5,  
                        A = 0,  
                        T = 262144
                        )
  tName <- names(transforms)[[i]]
  tName <- gsub(" ", "_", tName)
  tName <- paste0("tr", tName)
  t@parameters <- compensatedParameter(parameters = names(transforms)[[i]],
                            spillRefId = "compMat",
                            transformationId = tName,
                            searchEnv = flowEnv)
  t@transformationId <- tName
  
  transforms[[i]] <- t
  names(transforms)[[i]] <- tName
  
  assign(names(transforms)[[i]], transforms[[i]])
}

# save our compensation matrix to the env
assign("compMat", compMat, envir=flowEnv)


# MAIN GATE PROCESSING LOOP
# loop through gates, attach compensation and transformation info
# flowUtils will automagically create transformations from compensation matrices
# we rescale gate boundaries so flowJo can display them properly

un <- 1;
for (i in 1:length(gates)){
  gateName <- names(gates)[[i]]
  message("\n", gateName)
  # rescale gate boundaries since flowJo applies gates on untransformed data 

  if (is(gates[[i]], "polygonGate")){
     message("polygonGate ")
     mm <- gates[[i]]@boundaries
     for (x in 1:ncol(mm)){
       if (colnames(mm)[x]  %in% scaleme){ # this gate is not FSC/SSC so must be scaled
          message(colnames(mm)[x], " rescaling for flowJo")
          # fj <- flowJoTrans(channelRange = max(mm[,x]), inverse = T)
          fj <- flowJoTrans(channelRange = markers[colnames(mm)[x],]$range, 
                            #maxValue = markers[colnames(mm)[i],]$maxRange, 
                            inverse = T)
          
          mm[,x] <- fj(mm[,x])
       }
     }
     gates[[i]]@boundaries <- mm
    # writeGatingML doesn't use channel names in PolygonGate@boundaries 
    colnames(gates[[i]]@boundaries) <- NULL
    
  } else if (is(gates[[i]], "rectangleGate") ){ 
    # NOTE: rectangleGate can have -Inf/Inf coordinates
    # so check and only transform finite values with flowJoTrans
     message("rectangleGate ")
     mm <- gates[[i]]@min
     for (x in 1:length(gates[[i]]@min)){
       if (names(mm)[x] %in% scaleme && is.finite(mm[x]) ){
         message(names(mm)[x], " rescaling min for flowJo")
         # fj <- flowJoTrans(channelRange = max(mm[x]), inverse = T)
         fj <- flowJoTrans(channelRange = markers[names(mm)[x],]$range, 
                           #maxValue = markers[names(mm)[x],]$maxRange, 
                           inverse = T)
         
         mm[x] <- fj(mm[x])
       }     
     }
    gates[[i]]@min <- mm
     
    mm <- gates[[i]]@max
     for (x in 1:length(gates[[i]]@min)){
        if (names(mm)[x] %in% scaleme && is.finite(mm[x]) ){
          message(names(mm)[x], " rescaling max for flowJo")
          fj <- flowJoTrans(channelRange = max(mm[x]), inverse = T)
          mm[x] <- fj(mm[x])
        }

  }
    gates[[i]]@max <- mm
    
  } else if (is(gates[[i]], "ellipsoidGate")){
    message("ellipsoidGate ")
    
  }
  
  len <- length(parameters(gates[[i]]))
  trPars <- vector("list", len)
  for (j in 1:len){ # reference the transformations in gate parameters
    chanl <- eval(gates[[i]]@parameters[[j]]@parameters)
    if (chanl %in% rownames(markers) && !is.na(markers[chanl,"desc"])){ # have compensation, use it
      #trPars[j] <- unlist(transforms[which(names(transforms) == chanl)])
      tr <- gsub(" ", "_", chanl)
      tr <- paste0("tr", tr)
      tr <- get(tr)
      tr@parameters <- compensatedParameter(parameters = chanl,
                                            spillRefId = "compMat",
                                            transformationId = chanl,
                                            searchEnv = flowEnv)
      trPars[[j]] <- tr 
      message(chanl, " compensated")

    } 
   
    else { # uncompensated gate, define a unitytransform which does nothing to data
      trPars[[j]] <- unitytransform(chanl, transformationId = paste0("unity",un))
      un <- un + 1
      message(chanl, " no compensation!")
    }
    
  }
  gates[[i]]@parameters <- new("parameters", trPars)  

} # end main gate processing loop


# create parent <-> child gate relationships
# names(parents)[i] is gate, parents[i] is that gate's parent

for (i in 1:length(parents)){
  gateName <- gateNames[[i]]
  parent <- parents[[i]]
  if (parent == "root"){ # top-most gate, just write the gate
      child <- gates[[i]]
      gates[[i]] <- child
      flowEnv[[gateName]] <- gates[[i]]
      
  } else { # use subsetFilter to refer to parent
    child <- new("subsetFilter", 
                 filters=list(gates[[i]], flowEnv[[parent]]), 
                 filterId=gateName
          )
    gates[[i]] <- child
    flowEnv[[gateName]] <- gates[[i]]
    # strip the filterId of the leaf gate so we don't write extraneous XML
    flowEnv[[gateName]]@filters[[1]]@filterId <- ""
  }
  
  
  # go through the environment and strip out the children's filterID 
  
}

# NO HIERARCHY!!!
# for (i in 1:length(gates)){
#   gateName <- gateNames[[i]]
#   flowEnv[[gateName]] <- gates[[i]]
# }

#child <- new("subsetFilter", filters=list(child, flowEnv[[parent]]) )


# write the environment out to file
write.gatingML(flowEnv, outFile)

# for debugging
# flowEnvOut <<- flowEnv

message("\n Successfully wrote GatingML file ", outFile)

} # end function definition
