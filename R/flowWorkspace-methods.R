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

##################################################
#merge GatingSets into groups based on their gating schemes
#Be careful that the splitted resluts still points to the original data set!!
#drop the unused channels if needed before merging them 
##########################################################

#merge GatingSetList to GatingSet

#from list to GatingSetList
setMethod("merge",signature=c("list"),function(x,...){
#      browser()
      
      message("Grouping by Gating tree...")
      node_seq <-unlist(lapply(x,function(this_gs){
                this_gh <- this_gs[[1]]
                this_nodes <- getNodes(this_gh,isPath=T)
                paste(this_nodes,collapse = "")
                
              }))
      gs_groups <- split(x,node_seq)
        
      
      #start to merge      
      lapply(1:length(gs_groups),function(i){
#            browser()
            this_group <- gs_groups[[i]]

            #drop the unused marker from fs
            if(length(this_group) > 1){
              this_group <- lapply(this_group,function(this_gs){
#                    browser()
                    this_fs <- getData(this_gs)
                    
                    this_pd <- pData(parameters(this_fs[[1]]))
                    non_na_channel <- unname(!is.na(this_pd[,"desc"]))
                    to_include <- grepl(pattern="[FS]SC|[Tt]ime",this_pd[,"name"])
                    to_include <- to_include |  non_na_channel
                    if(length(which(to_include)) != nrow(this_pd)){
                    
                      message("drop empty channel:",this_pd[!to_include,1])
                      
                      ncFlowSet(this_gs) <- this_fs[,to_include]
                    
                    }
                    this_gs
                  })
            }
            
            GatingSetList(this_group)
          })
        
      
    })
      
######################################
##archive/unarchive to/from a folder 
##it is faster than tar-version,but require
##a new dest folder to avoid overwriting
##the old data by mistake
##currently not exposed to end user
######################################
save_gs<-function(G,path,overwrite = FALSE,...){
#  browser()
  
  if(file.exists(path)){
    path <- normalizePath(path,mustWork = TRUE)
    if(overwrite){
#      browser()
      res <- unlink(path, recursive = TRUE)
      if(res == 1){
        stop("failed to delete ",path)
      }
    }else{
      stop(path,"' already exists!")  
    }
    
  }
  
  dir.create(path = path)
  #do the dir normalization again after it is created
  path <- normalizePath(path,mustWork = TRUE)
  invisible(flowWorkspace:::.save_gs(G,path = path, ...))
  message("Done\nTo reload it, use 'load_gs' function\n")
  
  
}


load_gs<-function(path){
#  browser()
  path <- normalizePath(path,mustWork = TRUE)
  if(!file.exists(path))
    stop(path,"' not found!")
  files<-list.files(path)
#   browser()
  flowWorkspace:::.load_gs(output = path, files = files)$gs
  
}