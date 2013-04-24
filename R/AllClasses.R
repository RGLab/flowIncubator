
setClass("GatingSetList"
        ,representation=representation(
            data = "list"
            ,samples="character" #this determine the order of samples exposed to user
          )
          ,validity=function(object){
            
            gs_list <- object@data
            #check overlapping samples
            gs_samples <- unlist(lapply(gs_list,getSamples))
            if(any(duplicated(gs_samples))){
              return ("There are overlapping samples across GatingSets!")
            }
            
            
            gs1 <- gs_list[[1]]
            
            #compare GatingSets
            
            res <- sapply(gs_list[-1],function(this_gs){
                  
                  .compareGatingSet(this_gs,gs1)
                })
            
            
            is_error <- sapply(res,function(this_res){
                            class(this_res) == "character"
                          })
#            browser()
            if(any(is_error)){
              this_error_ind <- which(is_error)[1]
              return (paste("GatingSet 1 and",this_error_ind+1,":",res[this_error_ind]))
            }
            #check sample vector
            if(!.isValidSamples(object@samples,gs_list)){
              return ("'samples' slot is not consisitent with sample names from GatingSets!")
            }          
            return (TRUE)
          }
        )
.flattenedGatingHiearchy<-function(gh){
  this_nodes <- getNodes(gh,isPath=T)
  paste(this_nodes,collapse = "")
}        
#TODO:gating tree comparison needs to be improved        
.compareGatingHierarchy<-function(gh1,gh2){
  if(identical(.flattenedGatingHiearchy(gh1),.flattenedGatingHiearchy(gh2))){
    return (TRUE)
  }else{
    return (paste("gating structure doesn't match:",getSample(gh1),getSample(gh2)))
  }
}
.compareFlowData<-function(fs1,fs2){
  col1 <- colnames(fs1)
  col2 <- colnames(fs2)
  if(!identical(col1,col2)){
#    sCol1 <- paste(col1,collapse=" ")
#    sCol2 <- paste(col2,collapse=" ")
    msg <- paste("colnames of flowSets don't match!")
    return (msg)
  }
  if(!identical(colnames(pData(fs1)),colnames(pData(fs2)))){
    return ("pData of flow data doesn't match!")
  }
  return  (TRUE)
  
}
.compareGatingSet<-function(gs1,gs2){
  gh1 <- gs1[[1]]
  gh2 <- gs2[[1]]
  res <- .compareGatingHierarchy(gh1,gh2)
  if(class(res) == "character"){
    return (res)
  }
  fs1 <- getData(gs1)
  fs2 <- getData(gs2)
  .compareFlowData(fs1,fs2)
}
#validity check for samples slot        
.isValidSamples<-function(samples,object){
  
  return (setequal(unlist(lapply(object,getSamples)),samples))
}
          
## Constructor
GatingSetList <- function(x,samples = NULL)
{
  names(x)<-NULL#strip names from the list because rbind2 doesn't like it
  flowCore:::checkClass(x, "list")
  if(is.null(samples)){
    samples <- unlist(lapply(x,getSamples))
  }
  x <- new("GatingSetList", data = x, samples = samples)
  return(x)
}