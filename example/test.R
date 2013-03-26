# TODO: Add comment
# 
# Author: wjiang2
###############################################################################


library(flowWorkspace)

gs <- flowWorkspace:::load_gs("/home/wjiang2/rglab/workspace/analysis/HVTN080/output/HVTN_small/")

getData(gs[[1]],"4+")
getData(gs[[1]],6)

getData(gs[1])
res <- getData(gs[1],quote(`4+/TNFa+|4+/IL2+`))
res[[1]]

getNodes(gs[[1]])
getData(gs[1],"4+")[[1]]
getData(gs[1],6)[[1]]

getData(gs[[1]])[getIndices(gs[[1]],"TNFa+")|getIndices(gs[[1]],"IL2+"),]
