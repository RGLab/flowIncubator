# TODO: Add comment
# 
# Author: wjiang2
###############################################################################
#library(flowWorkspace)
library(flowIncubator)


#new_stats <-unlist(sapply(getNodes(gs[[1]]),function(this_node)length(which(getIndices(gs[[1]],this_node)))))[1:20]
#cbind(getPopStats(gs[[1]])[1:20,3,drop=F],new_stats)

#getData
gs <- load_gs("~/rglab/workspace/analysis/HVTN080/output/HVTNsubset")

save_gs(gs,"~/rglab/workspace/flowIncubator/output/test")
gs <- load_gs("~/rglab/workspace/flowIncubator/output/test")
getData(gs)
save_gs(gs,path="~/rglab/workspace/flowIncubator/output/remote/test2",overwrite=T,cdf="symlink")
gs <- load_gs("~/rglab/workspace/flowIncubator/output/remote/test2")
save_gs(gs,path="~/rglab/workspace/flowIncubator/output/test d", cdf="move")
gh <- gs[[1]]
x11()
plotGate(gs[[1]],"4+/IL2+")
g <- getGate(gs[[1]],"4+/IL2+")
g@boundaries[c(1,4),1] <- c(2.8,2.8)
setGate(gs[[1]],"4+/IL2+",g)
getData(gs[[1]],"4+")
getData(gs[[1]],6)

getNodes(gh)

getNodes(gh,147)
setNode(gh,147,"d/test")
getProp(gh,"4+") 
getProp(gh,"TNFa+") #refer to 4+/TNFa+
getProp(gh,"144.TNFa+") #refer to 8+/TNFa+

getProp(gh,"3+/4+") 
getProp(gh,"4+/TNFa+")
getProp(gh,"8+/TNFa+") #TODO:find the first match of starting node

getProp(gh,"4+/IFNg+")
getProp(gh,"Not 4+/IFNg+")

Rm("3+/4+",gs)
getData(gs[1])
res_ind <- getIndices(gs[1],quote(`4+/TNFa+|4+/IL2+`))
indMat <- getIndiceMat(gs[[1]],quote(`4+/TNFa+|4+/IL2+|/4+/IFNg+`))
res <- getData(gs[1],quote(`4+/TNFa+|4+/IL2+|/4+/IFNg+`))



plot(gh)
getNodes(gh)
res[[1]]

getNodes(gs[[1]])
getData(gs[1],"4+")[[1]]
getData(gs[1],6)[[1]]

getData(gs[[1]])[getIndices(gs[[1]],"TNFa+")|getIndices(gs[[1]],"IL2+"),]


showMethods("[[",classes="GatingSetList")

library(flowIncubator)
#unloadNamespace("flowIncubator")

#load several GatingSets from disk
gs_list<-lapply(list.dirs("~/rglab/workspace/flowIncubator/output/gs_toMerge",full=T,recur=F)
              ,function(this_folder){
      load_gs(this_folder)
    })
plot(gs_list[[1]][[1]])
gslist2 <- GatingSetList(gs_list[c(1,3)])
pData(ncFlowSet(gs_list[[2]]))$id=NULL
pData(ncFlowSet(gs_list[[5]]))

gh<-gs_list[[1]][[1]]
getPopStats(gh)[,2:3]
setNode(gs_list[[1]],6,"time")
this_g <- getGate(gh,6)
str(this_g)
this_g@boundaries[,1] <- this_g@boundaries[,1] *100 
setGate(gh,6,this_g,negated=T)
recompute(gs_list[[1]])
x11()
plotGate(gh,6,xbin=64,xlim=c(-1000,10000))
xyplot(`<Pacific Blue-A>`~`Time`,getData(gh),filter=getGate(gh,6),smooth=F,xlim=c(-1000,10000))
#gs_list is a list
gs_groups <- merge_gs(gs_list)
#returns a list of GatingSetList objects
gslist2 <- gs_groups[[1]]
#gslist2 is a GatingSetList that contains multiple GatingSets and they share the same gating and data structure
gslist2
class(gslist2)
getSamples(gslist2)


#reference a GatingSet by numeric index
gslist2[[1]]
#reference a GatingSet by character index
gslist2[["30104.fcs"]]

#loop through all GatingSets within GatingSetList
lapply(gslist2,getSamples)

#subset a GatingSetList by [
getSamples(gslist2[c(4,1)])
getSamples(gslist2[c(1,4)])
gslist2[c("30104.fcs")]

#get flow data from it
getData(gslist2)
#get gated flow data from a particular popoulation (by numeric or character index)
getData(gslist2,4)


#extract the gates associated with one popoulation
getGate(gslist2,"3+")
getGate(gslist2,5)


#extract the pheno data
pData(gslist2[2:1])
#modify the pheno data
pd <- pData(gslist2)
pd$id <- 1:nrow(pd)
pData(gslist2) <- pd
pData(gslist2[3:2])

#plot the gate
plotGate(gslist2[1:2],5,smooth=T)
plotGate_labkey(gslist2[3:4],4,x="<APC Cy7-A>",y="<PE Tx RD-A>",smooth=T)

#remove cerntain gates by loop through GatingSets
getNodes(gslist2[[1]])
lapply(gslist2,function(gs)Rm("Excl",gs))

#extract the stats
getPopStats(gslist2)
#extract statistics by using getQAStats defined in QUALIFIER package
res<-getQAStats(gslist2[c(4,2)],isMFI=F,isSpike=F,nslaves=1)


#archive the GatingSetList
save_gslist(gslist2, path ="~/rglab/workspace/flowIncubator/output/gslist",overwrite=T)
gslist2 <- load_gslist(path ="~/rglab/workspace/flowIncubator/output/gslist")

#convert GatingSetList into one GatingSet by rbind2
gs_merged2 <- rbind2(gslist2,ncdfFile=path.expand(tempfile(tmpdir="~/rglab/workspace/flowIncubator/output/",fileext=".nc")))
gs_merged2

gh<-gslist2[[1]]
getNodes(gh)
plot(gh,5)
plot(gh,"3+")
plot(gh,"19+ 20+")
plot(gh,"19+ 20-")
plot(gh,"19- 20-")
dev.off()
g<-flowWorkspace:::.getGraph()
plot(subGraph(nodes(g)[3:4],g))
nodes(g)
plot(g)


gh<-gs_list[[1]][[1]]

g<-getGate(gh,"Excl")
g@boundaries[,1]<-g@boundaries[,1]*100

str(getGate(gh,"Excl"))

xyplot(`<Pacific Blue-A>`~Time,getData(gh,"3+")
,filter=g
,xlim=c(-1000,10000)
,smooth=F
,overlay=getData(gh,"IFNg+")
)
getNodes(gh)
range(exprs(getData(gh,"3+"))[,1])


load("/loc/no-backup/ramey/Lyoplate3-fs_list.RData")
list1 <- lapply(fs_list,GatingSet)
gslist <- GatingSetList(list1)
gs <- rbind2(gslist1)


#merge by dropping redundant terminal gates
gs <- load_gs(file.path("/home/wjiang2/rglab/workspace/analysis/ITN507ST/Newell/sharedData","autoGating12"))
gs1<-clone(gs[1:2])
gs2<-clone(gs[3:4])
setNode(gs2,"CD3","cd3")
getNodes(gs2[[1]])
Rm("IgD+cd27+",gs2)
gs_groups <- .groupByTree(list(gs1,gs2))
toRemove <- .checkRedundantNodes(gs_groups)
.dropRedundantNodes(gs_groups,toRemove)
new_group <- unlist(gs_groups)
gslist <- .mergeGS(new_group)
