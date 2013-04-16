
setMethod("rbind2",
    signature=signature("GatingSetList","missing"),
    definition=function(x,y="missing",...)
    {
#           browser()
      isNcdfList<-lapply(x,function(gs)flowWorkspace:::isNcdf(gs[[1]]))
      if(all(duplicated(unlist(isNcdfList))[-1])){
#               browser()
        #combine flowset/ncdfFlowSet
        fsList<-lapply(x,getData)
        if(isNcdfList[[1]])
          fs<-rbind2(as(fsList,"ncdfFlowList"),...)
        else
        {
          ##using original flowCore::rbind2 for flowSet
          fs<-fsList[[1]]
          for(i in 2:length(fsList))
            fs<-rbind2(fs,fsList[[i]])
        }
        
        #combine tree structure
        ptrlist<-lapply(x,function(gs)gs@pointer)
        sampleList<-lapply(x,getSamples)
        pointer<-.Call("R_combineGatingSet",ptrlist,sampleList)
        G<-new("GatingSetInternal")
        G@pointer<-pointer
        
        #combine R objects
        ne<-new.env(parent=emptyenv());
        assign("ncfs",fs,envir=ne)
        set<-unlist(lapply(x,function(gs)gs@set))
        #deep copying of tree
        for(i in seq_along(set))
        {
          #create new local data environment that stores axis and flowData environment
          localDataEnvOld<-nodeDataDefaults(set[[i]]@tree,"data")
          localDataEnv<-new.env(parent=emptyenv())
          copyEnv(localDataEnvOld,localDataEnv)
          #update flowData environment with new ncfs
          assign("data",ne,localDataEnv)
          #sync back to tree
          nodeDataDefaults(set[[i]]@tree,"data")<-localDataEnv
          #upodate pointer
          set[[i]]@pointer<-pointer
        }
        
        G@set<-set
        
      }else{
        stop("Can't combine gating sets. They should all use the same storage method. (Netcdf, or not..)")
      }
      return(G);  
      
    })
setMethod("show",
    signature = signature(object="GatingSetList"),
    definition = function(object) { 
      cat("A GatingSetList with", length(object@data),"GatingSet\n")
      cat("containing", length(unique(getSamples(object))), " unique samples.") 
      cat("\n")
    })
setMethod("getSamples", 
    signature = signature(x = "GatingSetList"),
    function(x,...) {
      unlist(lapply(x@data,getSamples))      
    })
setMethod("lapply","GatingSetList",function(X,FUN,...){
      lapply(X@data,FUN,...)
    })
setMethod("[[",c(x="GatingSetList",i="numeric"),function(x,i,j,...){
      #convert non-character indices to character
#      browser()
      this_samples <- getSamples(x)
      nSamples <- length(this_samples)
      if(i > nSamples){
        stop(i, " is larger than the number of samples: ", nSamples)
      }
        x[[this_samples[i]]]
      
    })

setMethod("[[",c(x="GatingSetList",i="logical"),function(x,i,j,...){
      #convert non-character indices to character
      
      x[[getSamples(x)[i]]]
      
    })
setMethod("[[",c(x="GatingSetList",i="character"),function(x,i,j,...){
      #convert non-character indices to character
      gh <- NULL
      for(gs in x@data){
#              browser()
            this_samples <- getSamples(gs)
            ind <- match(i,this_samples)
            if(!is.na(ind)){
              gh <- gs[[ind]]
            }
      }
      if(is.null(gh)){
        stop(i, " not found in GatingSetList!")
      }else{
        return (gh)
      }
    })
setMethod("[",c(x="GatingSetList",i="numeric"),function(x,i,j,...){
#      browser()
      x[getSamples(x)[i]]
      
    })

setMethod("[",c(x="GatingSetList",i="logical"),function(x,i,j,...){
      
      x[getSamples(x)[i]]
   
    })
#TODO:add metaData slot to GatingSetList to maintain the sample order
setMethod("[",c(x="GatingSetList",i="character"),function(x,i,j,...){
#      browser()
      samples <- getSamples(x)
      noFound <- is.na(match(i,samples))
      if(any(noFound)){
        stop(i(noFound), "not found in GatingSetList!")
      }
      res <- lapply(x,function(gs){
#            browser()
                  this_samples <- getSamples(gs)
                  ind <- match(i,this_samples)
                  this_subset <- i[!is.na(ind)] 
                  if(length(this_subset)>0){
                    return (gs[this_subset])
                  }else{
                    NULL
                  }
                })
      res <- res[!unlist(lapply(res,is.null))]
      GatingSetList(res)
      
    })

setMethod("getData",c(obj="GatingSetList",y="missing"),function(obj,y,...){
      stop("node index 'y' is missing!")
    })
#
setMethod("getData",signature(obj="GatingSetList",y="numeric"),function(obj,y,...){
#      browser()
      this_node <- getNodes(obj[[1]])[y]
      getData(obj,this_node)
    })
setMethod("getData",c(obj="GatingSetList",y="character"),function(obj, y, max=30, ...){
#      browser()
      if(length(getSamples(obj))>max){
        stop("You are trying to return a flowSet for more than ", max, " samples!Try to increase this limit by specifing 'max' option if you have enough memory.")
      }
#      browser()
      res <- lapply(obj,function(gs){
                NcdfFlowSetToFlowSet(getData(gs,y))
          })
      fs<-res[[1]]
      for(i in 2:length(res))
        fs<-rbind2(fs,res[[i]])
      fs
    })

setMethod("pData","GatingSetList",function(object){
#      browser()
      res <- lapply(object,function(gs){
              pData(gs)      
          })
      do.call(rbind,res)
    })
setMethod("getGate",signature(obj="GatingSetList",y="numeric"),function(obj,y,tsort=FALSE){
      getGate(obj,getNodes(obj[[1]])[y])
    })
setMethod("getGate",signature(obj="GatingSetList",y="character"),function(obj,y,tsort=FALSE){
      res <- lapply(obj,function(gs){
            getGate(gs,y)      
          })
      unlist(res,recur=FALSE)
      
    })


setMethod("plotGate",signature(x="GatingSetList",y="numeric"),function(x,y, ...){
      selectMethod("plotGate",sig=c(x="GatingSetInternal",y="numeric"))(x=x, y=y, ...)
      
      
    })

setMethod("getQAStats",signature("GatingSetList"),function(obj,...){
      res <- lapply(obj,function(gs){
            getQAStats(gs)      
          })
      do.call(rbind2,res)
    })