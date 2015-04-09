#' Fetch a list of keywords from a GatingSet
#' Return them as a data.table with columns names for the keyword
#' The FIL keyword is renamed to 'name' for sample name consistency
#' No error checking at the moment.
#' @export 
setMethod("getKeywords",c("GatingSet","character"),function(obj,y){
      gs <- obj
      kv <- y

  r<-as.data.frame(do.call(cbind,lapply(kv,function(k){
    keyword(gs,k)[1]
  })))
  data.table::setnames(r,"$FIL","name")
  r
})
