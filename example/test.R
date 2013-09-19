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

gh<-gs_list[[1]][[1]]
getPopStats(gh)[,2:3]

x11()
plotGate(gh,6,xbin=64,xlim=c(-1000,10000))

#gs_list is a list
gslist <- flowIncubator:::.mergeGS(gs_list)

gslist
keyword(gslist@data[[1]],"Stim")
keyword(gslist[[1]],"Stim")
class(gslist)
getSamples(gslist)


#reference a GatingSet by numeric index
gslist[[1]]
#reference a GatingSet by character index
gslist2[["977457.fcs"]]

#loop through all GatingSets within GatingSetList
lapply(gslist,getSamples,level=1)

#subset a GatingSetList by [
getSamples(gslist[c(4,1)])
getSamples(gslist[c(1,4)])
gslist[c("977457.fcs")]

#get flow data from it
getData(gslist)
#get gated flow data from a particular popoulation (by numeric or character index)
getData(gslist,4)


#extract the gates associated with one popoulation
getGate(gslist,"3+")
getGate(gslist,5)


#extract the pheno data
pData(gslist[2:1])
#modify the pheno data
pd <- pData(gslist)
pd$fileid <- 1:nrow(pd)
pData(gslist) <- pd
pData(gslist[3:2])

#plot the gate
plotGate(gslist[1:2],5,smooth=T)
plotGate_labkey(gslist[3:4],4,x="<APC Cy7-A>",y="<PE Tx RD-A>",smooth=T)

#remove cerntain gates by loop through GatingSets
getNodes(gslist[[1]])
#lapply(gslist,function(gs)Rm("Excl",gs))

#extract the stats
getPopStats(gslist)
#extract statistics by using getQAStats defined in QUALIFIER package
res <- getQAStats(gslist[c(4,2)],isMFI=F,isSpike=F,nslaves=1)


#archive the GatingSetList
save_gslist(gslist, path ="~/rglab/workspace/flowIncubator/output/gslist",overwrite=T)
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
g<-flowWorkspace:::.getGraph(gh)
plot(graph::subGraph(nodes(g)[3:4],g))
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

library(flowIncubator)
#merge by dropping redundant terminal gates
gs <- load_gs(file.path("/home/wjiang2/rglab/workspace/analysis/ITN507ST/Newell/sharedData","autoGating12"))
gs1<-clone(gs[1:2])
gs2<-clone(gs[3:4])
setNode(gs2,"CD3","cd3")
getNodes(gs2[[1]])
Rm("IgD+cd27+",gs2)


gs_groups <- flowIncubator:::.groupByTree(list(gs1,gs2))
toRemove <- flowIncubator:::.checkRedundantNodes(gs_groups)
flowIncubator:::.dropRedundantNodes(gs_groups,toRemove)
new_group <- unlist(gs_groups)
gslist <- flowIncubator:::.mergeGS(new_group)



#impute Gate
getNodes(gs[[1]])
Rm("MTG_gate", gs)
library(Cairo)
CairoX11()
#visual check gating
plotGate(gs, "MTG_gate", xbin=64,margin=T
    ,type="densityplot"
    ,darg=list(bw = "nrd0",n=512)
)
#outlier detection
popStats <- getPopStats(gs)["/s1/s2/live/lymph/cd3/Bcell/MTG_gate",]
outlierInd <- QUALIFIER:::outlier.cutoff(popStats,uBound=0.80)
failedSamples <- names(which(outlierInd))

hist(popStats)
#.nearestSample(gs, "MTG_gate", target = "M+T panel_903997-25.fcs", source = getSamples(gs)[-5])
system.time(refSamples <- .nearestSamples(gs, "MTG_gate", failedSamples, method = "ks.test"))
samplePairs <- data.frame(fail=names(refSamples),ref=refSamples,row.names=NULL)
refGates <- sapply(refSamples,function(i)getGate(gs[[i]],"MTG_gate"))
setGate(gs[failedSamples],"MTG_gate",refGates)
recompute(gs[failedSamples],"MTG_gate")

gslist <- load_gslist("/loc/no-backup/Leo/338ac112e41bc368f3d963e31a374f04/")
gslist <- save_gslist_labkey(gslist, path = "/loc/no-backup/mike/tmp3",overwrite=T, cdf = "move")
gslist <- load_gslist("/loc/no-backup/mike/tmp3")
getData(gslist[[1]])


