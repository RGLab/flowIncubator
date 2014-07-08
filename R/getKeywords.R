#' Fetch a list of keywords from a GatingSet
#' Return them as a data.table with columns names for the keyword
#' The FIL keyword is renamed to 'name' for sample name consistency
#' No error checking at the moment.
getKeywords<-function(gs,kv){
  r<-as.data.frame(do.call(cbind,lapply(kv,function(k){
    keyword(gs,k)[1]
  })))
  data.table::setnames(r,"$FIL","name")
  r
}
