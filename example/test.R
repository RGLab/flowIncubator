# TODO: Add comment
# 
# Author: wjiang2
###############################################################################
library(flowWorkspace)
library(flowIncubator)


new_stats <-unlist(sapply(getNodes(gs[[1]]),function(this_node)length(which(getIndices(gs[[1]],this_node)))))[1:20]
cbind(getPopStats(gs[[1]])[1:20,3,drop=F],new_stats)

#getData
gs <- load_gs("~/rglab/workspace/analysis/HVTN080/output/HVTNsubset1")
gs <- load_gs("/home/wjiang2/rglab/workspace/analysis/HVTN080/output/HVTN_small/")
x11()
plotGate(gs[[1]],smooth=TRUE)
getData(gs[[1]],"4+")
getData(gs[[1]],6)

getData(gs[1])
res <- getData(gs[1],quote(`4+/TNFa+|4+/IL2+`))
res[[1]]

getNodes(gs[[1]])
getData(gs[1],"4+")[[1]]
getData(gs[1],6)[[1]]

getData(gs[[1]])[getIndices(gs[[1]],"TNFa+")|getIndices(gs[[1]],"IL2+"),]



#split and merge
gs_list<-lapply(list.files("~/rglab/workspace/flowIncubator/output/gs_toMerge",full=T),function(this_folder){
      load_gs(this_folder)
    })


gs_groups <- split(GatingSetList(gs_list))

gs_list_merged <- merge(GatingSetList(gs_list),path="flowIncubator/output/")
getwd()
gs_list_merged

getData(gs_list_merged[[1]],quote(`4+/TNFa+|4+/IL2+`))[[1]]
getData(gs_list_merged[[1]])

