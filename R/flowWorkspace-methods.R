#routine to return the indices by specify boolean combination of reference nodes:
# y is a quoted expression.
#1.adds the boolean gates and does the gating on the fly 
#2.return the indices associated with that bool gate
#3. remove the bool gate
# the typical use case would be extracting any-cytokine-expressed cells
#example:
# getIndices(gs,quote(`4+/TNFa+|4+/IL2+`)) (it may be faster than R version
###############################################################################
setMethod("getIndices",signature=c("GatingSetInternal","name"),function(obj, y, ...){
      
      bf <- eval(substitute(booleanFilter(v),list(v=y)))
      gh <- obj[[1]]
      
          suppressMessages({
              suppressWarnings(
                id <- add(obj,bf)
                )
                
              allNodes <- getNodes(gh,isPath=TRUE)
              this_node <- allNodes[id]
              
              
              res <-try(recompute(obj,id),silent=T)
           })
      
      
      if(class(res)=="try-error"){
        Rm(this_node,obj)
        stop(res)
      }else{
        this_ind <- lapply(obj,function(this_gh)getIndices(this_gh,this_node))
        Rm(this_node,obj)
        this_ind  
      }
      
   })

getIndiceMat<-function(gh,y){
  strExpr <- as.character(y)
  nodes <- strsplit(strExpr,split="\\|")[[1]]
#  browser()
  #extract logical indices for each cytokine gate
  indice_list <- sapply(nodes,function(this_node)getIndices(gh,this_node),simplify = FALSE)
  #construct the indice matrix
  do.call(cbind,indice_list)
}

setMethod("getData",signature=c("GatingSetInternal","name"),function(obj, y,...){
      #get ind of bool gate
      bool_inds <- getIndices(obj,y,...)
      
      #########################################
      #create mapping between pops and channels
      ##########################################
      #get pop names
      strExpr <- as.character(y)
      popNames <- strsplit(strExpr,split="\\|")[[1]]
      
      #parse the markers of interest from pop names
      markers_selected <- sapply(popNames,function(this_pop){
            this_pops <- strsplit(split="/",this_pop)[[1]]
            #get the terminal node
            term_pop <- this_pops[length(this_pops)]
            term_pop
          },USE.NAMES=FALSE)
      
      #match to the pdata of flow frame
      fr <- getData(obj[[1]])
      this_pd <- pData(parameters(fr))
      all_markers <- this_pd[,"desc"]
      all_markers <- as.character(all_markers)
      all_markers[is.na(all_markers)] <- "NA"
      is_matched <- sapply(all_markers,function(this_marker){
            this_matched <- grep(pattern = this_marker, x=markers_selected, fixed=TRUE)
            if(length(this_matched)>1){
              stop("multiple populations mached to:", this_marker)
            }else if(length(this_matched)==0){
              res <- FALSE
            }else{
              res <- TRUE
              names(res) <- popNames[this_matched]
            }
            res
          }, USE.NAMES=FALSE)
#      browser()
      #get mapping
      pop_chnl <- cbind(pop=names(is_matched[is_matched]),this_pd[is_matched,c("name","desc")])
      
      
      lapply(obj,function(gh){
            #get mask mat
#      browser()
            this_sample <- getSample(gh)
            this_ind <-  bool_inds[[this_sample]]
            this_pops <-  pop_chnl[,"pop"]
            this_mat <- getIndiceMat(gh,y)[this_ind,this_pops]
            #subset data by channels selected
            this_chnls <- pop_chnl[,"name"]
            this_data <- getData(gh)
            this_subset <- exprs(this_data)[this_ind,this_chnls] 
            #masking the data
            this_subset <- this_subset *  this_mat
            colnames(this_subset) <- pop_chnl[,"desc"]
            this_subset
          })
      
    })
      
##################################################
#merge GatingSets into groups based on their gating schemes
#Be careful that the splitted resluts still points to the original data set!!
#drop the unused channels if needed before merging them 
##########################################################


#from list to GatingSetList
merge_gs<-function(x,...){
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
        
      
    }
      
######################################
##archive/unarchive to/from a folder 
##it is faster than tar-version,but require
##a new dest folder to avoid overwriting
##the old data by mistake
##currently not exposed to end user
######################################
save_gs<-function(G,path,overwrite = FALSE, save.cdf = TRUE, ...){
#  browser()
  guid <- gs@guid
  rds_toSave <- paste(guid,"rds",sep=".")
  dat_toSave <- paste(guid,"dat",sep=".")
  
  if(file.exists(path)){
    path <- normalizePath(path,mustWork = TRUE)
    if(overwrite){
      this_files <- list.files(path)
      #validity check for non-empty folder
      if(length(this_files)!=0)
      {
        rds_ind <- grep("\\.rds$",this_files)
        dat_ind <- grep("\\.dat$",this_files)
        
        if(length(rds_ind)!=1||length(dat_ind)!=1){
          stop("Not a valid GatingSet archiving folder!")
        }else{
          this_rds <- this_files[rds_ind]
          this_dat <- this_files[dat_ind]
          
          if(this_rds!=rds_toSave||this_dat!=dat_toSave){
            stop("The GatingSet doesn't match the archived files in: ", path)
          }
        }
      }
      
      #validity check for cdf
      if(flowWorkspace:::isNcdf(G[[1]])){
        if(length(this_files)!=0){
          cdf_ind <- grep("\\.nc$",this_files)
          if(length(cdf_ind) != 1){
            stop("Not a valid GatingSet archiving folder!")
          }  
        }
        
      }
      if(length(this_files)!=0)
      {
        #start to delete the old files in path
        file.remove(file.path(path,rds_toSave))
        file.remove(file.path(path,dat_toSave))
        if(flowWorkspace:::isNcdf(G[[1]])&&save.cdf){
          file.remove(file.path(path,this_files[cdf_ind]))
        }
      }
      
    }else{
      stop(path,"' already exists!")  
    }
    
  }else{
    dir.create(path = path)
    #do the dir normalization again after it is created
    path <- normalizePath(path,mustWork = TRUE)
    
  }
#  browser()
  invisible(flowWorkspace:::.save_gs(G,path = path, copy.cdf = save.cdf, ...))
  message("Done\nTo reload it, use 'load_gs' function\n")
  
  
}


load_gs<-function(path){
#  browser()
  path <- normalizePath(path,mustWork = TRUE)
  if(!file.exists(path))
    stop(path,"' not found!")
  files<-list.files(path)
#   browser()
  gs <- flowWorkspace:::.load_gs(output = path, files = files)$gs
  guid <- try(slot(gs,"guid"),silent=T)
  if(class(guid)=="try-error"){
    #generate the guid for the old archive
    gs@guid <- system("uuidgen",intern = TRUE)
  }
  gs
}

.getAllDescendants <- function(gh,startNode,nodelist){
  
  children_nodes <- getChildren(gh,startNode)
  if(length(children_nodes)>0){
    for(this_parent in children_nodes){
      nodelist$v <- c(nodelist$v, this_parent)
      .getAllDescendants (gh,this_parent,nodelist)
    }  
  }
  
}
#plot subgraph
#TODO:merge with plot method in flowWorkspace
setMethod("plot",c("GatingHierarchyInternal","numeric"),function(x,y,...){
      
      
      # get graphNEL object
      g <- flowWorkspace:::.getGraph(x)
      nodelist <- new.env(parent=emptyenv())
      nodelist$v <-integer()
      .getAllDescendants (x,y,nodelist)
      nodelist$v <- c(nodelist$v,y)
      #assume the number y is consistent with  R graph node name: N_x 
      subNodes <- paste("N",nodelist$v-1,sep="_")
      #convert numeric index to node name
#      allNodes <- nodes(g)
#      startNode <- paste("N",y,sep="_")
     
      #get ride of parents (since dfs doesn't like it)
#      g <- subGraph(allNodes[-(1:(y-1))],g)
      #do dfs search
#      subGraphs <-dfs(g,startNode)
#      if(any(grepl("discovered",names(subGraphs)))){
#        
#        subNodes <- subGraphs$discovered  
#      }else{
#        #if muliple graphs,pick the first one that contains the current startNode
#        subGraphs <- subGraphs[[1]]
#        subNodes <- subGraphs$discovered
#      }
      
      if(length(subNodes)<=1){
        stop("Rgraphviz doesn't know how to plot leaf node!")
      }
      g <- subGraph(subNodes, g)
      
      flowWorkspace:::.plotGatingTree(g,...)

    })

setMethod("plot",c("GatingHierarchyInternal","character"),function(x,y,...){
      
      plot(x,flowWorkspace:::.getNodeInd(x,y))
      
    })