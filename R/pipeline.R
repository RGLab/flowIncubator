#' Step 1: determine which group to parse for each xml workspace
#' 
#' The regular expression may not always pick up the correct group for each xml.
#' Thus it is safer to separate this step from parsing so that users get to review and edit(if necessary) 
#' the resolved group results.
#'   
#' 
#' @param files character vector, full paths of `xml` workspace files
#' @param pattern character a regular expression to specifiy which group to parse
#' @return a list of opened workspace objects and the group names matched by regular expression.
#' @export 
#' @examples 
#' \dontrun{
#'  wsFiles <- list.files(dataPath, pattern = ".xml", full = T)
#'  groups <- get_groups(wsFiles, "(_MD)|(_SK)")
#' }
get_groups <- function(files, pattern){
  lapply(files, function(wsFile){
    
    
        batchID <- basename(wsFile)
        message(batchID)
        
        batchID <- sub(".xml", "", batchID)
        # parse workspace
        ws <- openWorkspace(wsFile, options = 1)
        
        sg <- getSampleGroups(ws)
        groupNames <- as.character(unique(sg$groupName)) 
        gn <- groupNames[grep(pattern, groupNames)]
        list(ws = ws , group = gn)
        
      })
    
}

#' Step 2: parse the entire batches
#' 
#' Sub-steps:
#' 1. parse workspaces into GatingSets 
#' 2. save the parsed GatingSets into "parsed" subfolder under output.dir 
#' 3. report the possible discrepancy within each GatingSet  
#' 
#' @param resolved_groups the results returned by 1st step
#' @param ... other arguments passed to 'parseWorkspace'
#' @return nothing
#' @export 
#' @examples 
#' \dontrun{
#'  
#'  
#' }
parse_xmls <- function(resolved_groups, output.dir, ...){
  
  
  
  lapply(resolved_groups, function(res){
        
        group <- res[["group"]]
        ws <- res[["ws"]]
        
        batchID <- basename(ws@file)
        batchID <- sub(".xml", "", batchID)
        message(batchID)
        
        sid <- unique(subset(sg, groupName == group)$sampleID)
        samps <- getSamples(ws)
        
        fcs <- unique(subset(samps, sampleID %in% sid)$name)
        
        
        gsName <- paste(batchID, "_" , gsub(' ','',group), sep = "")
        gsPath <- file.path(output.dir,"parsed",gsName)
        
        if(dir.exists(gsPath)){
          message("skip parsing ", batchID, "because ", gsPath, "already exists.")
        }else{
          
          
          ncdfFile <- file.path(output.dir,"output", paste(gsName,".nc",sep = ""))
          
          if(nrow(thisSamples)>0)
          {
            
            # parse
            res <- suppressMessages(try(gs <- parseWorkspace(ws
                                                , name = gn 
                                                , ncdfFile = ncdfFile
                                                , path = fcsPath
                                                , ...
                                            )))
            
            #check the discrepancy among gating trees (some time individual sample has extra terminal ndoe, which is mostly added by mistake)
            if(class(res) == "GatingSet"){
              #archive it
              save_gs(gs, path = gsPath, overwrite = TRUE, cdf = "move")
              gs_groups <- flowWorkspace:::.groupByTree(gs)
              return(length(gs_groups))
            }else
              return (res)
            
          }   
        }
  
  
  
      })
  
  
}

dropRedundantLeafNodes <- function(gs){
  gs_groups <- flowWorkspace:::.groupByTree(gs)
  gs_groups <- flowWorkspace:::.groupByTree(gs_groups)
  toRm <- flowWorkspace:::.checkRedundantNodes(gs_groups)
  flowWorkspace:::.dropRedundantNodes(gs_groups, toRm)
}

#' Step 3: merge the parsed gs into GatingSetList
#'
#' Sub-steps:
#' 1. load the GatingSets parsed from step.
#' 2 and does validity check to see if the gating trees within each gatingSet is clean
#'       try to resolve the discrepancy of gating trees within single GatingSet by dropping the extra terminal node for some samples (mistakenly added by the person who does flowJo)
#'      when there is non-leaf node causing the discrepancy, we throw the error to ask for human intervention such as: correct the typo of population name or drop the entire subtree
#' 3. drop the redundant channels (since they are usually one of the causes for channel discrepancy across batches)
#' 4. further drop the non-leaf nodes that cause the cross-batch discrepancy  
merge_gs <- function(output.dir){
  #' load parsed gs
  gs_list <- sapply(list.files(file.path(output.dir, "parsed"), recur = F, full = T),load_gs)
  gs_groups <- flowIncubator:::.groupByTree(gs_list)
}