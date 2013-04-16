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


showMethods("[",classes="GatingSetList")
library(flowWorkspace)
library(flowIncubator)
#merge
gs_list<-lapply(list.files("~/rglab/workspace/flowIncubator/output/gs_toMerge",full=T),function(this_folder){
      load_gs(this_folder)
    })


gs_groups <- merge(gs_list)

gslist2 <- gs_groups[[2]]
gslist2
class(gslist2)
getSamples(gslist2)

gslist2[[1]]
gslist2[["30104.fcs"]]
gslist2[1]
getSamples(gslist2[c(4,1)])
getSamples(gslist2[c(1,4)])
gslist2[c("30104.fcs")]
getData(gslist2)
getData(gslist2,4)
getGate(gslist2,"4+")
getGate(gslist2,6)
pData(gslist2[3:1])
plotGate(gslist2[1:2],7,smooth=T)
plotGate_labkey(gslist2[3:4],6,x="<APC Cy7-A>",y="<Alexa 680-A>",smooth=T)



pd <- pData(gslist2)
pd$id <- 1:nrow(pd)
pData(gslist2) <- pd
pData(gslist2[3:2])
getNodes(gslist2[[1]])
lapply(gslist2,function(gs)Rm("Excl",gs))
res<-getQAStats(gslist2[c(4,2)],isMFI=F,isSpike=F,nslaves=1)


gs_merged2 <- rbind2(gslist2,ncdfFile=path.expand(tempfile(tmpdir="~/rglab/workspace/flowIncubator/output/",fileext=".nc")))

gs_merged2

getData(gs_merged2[[1]],quote(`4+/TNFa+|4+/IL2+`))[[1]]
getData(gs_merged2[[1]])

