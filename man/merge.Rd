\name{merge}
\alias{merge}
\alias{merge-methods}
\alias{merge,GatingSetList-method}


\title{
merge GatingSets with identical gating trees.
}
\description{
merge GatingSets with identical gating trees.
}
\usage{
\S4method{merge}{GatingSetList}(x,path,...)

}
\arguments{
  \item{x}{
A \code{GatingSetList} object . 
}
  \item{path}{
\code{character} path to store the temporay cdf files 
}
}
\details{
group all the individual GatingSets based on their gating schemes by comparing the node list returned by \code{\link{getNodes}},which assumes
they follow the same population naming conventions. 
Meanwhile the unused channels are automatically dropped to make sure the \code{\link[flowWorkspace::rbind2]{rbind2}}
work successfully in the end.  
 }
\value{

A \code{\link{GatingSetList}} that contains multiple GatingSets each of which share the same gating and data structure.
}
\references{
\url{http://www.rglab.org/}
}
\author{
Mike Jiang \email{wjiang2@fhcrc.org}
}
\note{

}


\seealso{
\code{\link{rbind2}},\code{\link{GatingSetList}}
}

\examples{
\dontrun{
	#load gatingsets from disk
	#gs_toMerge is the path that stores multiple archived gatingsets
	gs_list<-lapply(list.files("flowIncubator/output/gs_toMerge",full=T),function(this_folder){
      flowWorkspace:::load_gs(this_folder)
    })

	#do the merging
	gs_list_merged <- merge(GatingSetList(gs_list))
	#print the merged result
	gs_list_merged


	
}
}
\keyword{ manip }
