#routine to return the data by specify boolean combination of reference nodes:
# y is a quoted expression.
#1.adds the boolean gates on the fly 
#2.return the data associated with that bool gate
#3. remove the bool gate
# the typical use case would be extracting any-cytokine-expressed cells
#example:
# getData(gs,quote(`4+/TNFa+|4+/IL2+`))
###############################################################################
setMethod("getData",signature=c("GatingSetInternal","name"),function(obj,y,...){
      
      bf <- eval(substitute(booleanFilter(v),list(v=y)))
      
      suppressMessages({
            id <- add(obj,bf)
            this_node <- getNodes(obj[[1]])[id]
            res <-try(recompute(obj,id),silent=T)
          })
      
      
      if(class(res)=="try-error"){
        Rm(this_node,obj)
        stop(res)
      }else{
        this_data <- getData(obj,this_node)
        Rm(this_node,obj)
        this_data  
      }
      
      
      
    })
