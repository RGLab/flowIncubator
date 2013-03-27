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

#split GatingSets into groups based on their gating schemes
#Be careful that the splitted resluts still points to the original data set!!
setMethod("split",signature=c("GatingSetList","missing"),function(x,f,...){
#      browser()
      node_seq <-unlist(lapply(x,function(this_gs){
                this_gh <- this_gs[[1]]
                this_nodes <- getNodes(this_gh,isPath=T)
                paste(this_nodes,collapse = "")
                
              }))
      split(x,node_seq)
      
    })
      


#drop the unused channels if needed before merging them 
setMethod("merge",signature=c("GatingSetList"),function(x,path = tempdir(),...){
#      browser()
      
      message("Grouping by Gating tree...")
      gs_groups <- split(x)
      res <- lapply(1:length(gs_groups),function(i){
#            browser()
            this_group <- gs_groups[[i]]
            if(length(this_group) == 1){
              return(clone(this_group[[1]]))
            }else{
              message("merging Group:",i)
              #drop the unused marker from fs
              this_group_new <- lapply(this_group,function(this_gs){
#                    browser()
                    this_fs <- getData(this_gs)
                    
                    this_pd <- pData(parameters(this_fs[[1]]))
                    non_na_channel <- unname(!is.na(this_pd[,"desc"]))
                    to_include <- grepl(pattern="[FS]SC|[Tt]ime",this_pd[,"name"])
                    to_include <- to_include |  non_na_channel
                    if(length(which(to_include)) != nrow(this_pd)){
                    
                      message("drop empty channel:",this_pd[!to_include,1])
                      this_folder <- path.expand(path)
                      #      browser()
                      this_fs <- clone.ncdfFlowSet(this_fs[,to_include]
                          ,isEmpty = FALSE
                          ,ncdfFile=tempfile(tmpdir=this_folder,fileext=".nc")
                      )
                      
                      #only update the data environment here to avoid changing original gs data by mistake
                      #since rbind2 is going to take care of the rest of environments
#                      browser()
                      gdata<-new.env(parent=emptyenv());
                      for(i in 1:length(this_gs@set)){
                        nd<-this_gs@set[[i]]@tree@nodeData
                        #local data env
                        nd@defaults$data<-new.env(hash=TRUE, parent=emptyenv())
                        #global data env
                        nd@defaults$data[["data"]]<-gdata
                        this_gs@set[[i]]@tree@nodeData<-nd
                      }
                      
                      ncFlowSet(this_gs) <- this_fs
                    }
                    this_gs
                  })
                 
              
              rbind2(GatingSetList(this_group_new))  
            }
            
          })
          
        message("merging finished!")
        GatingSetList(unlist(res))
    })
      
