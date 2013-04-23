\name{merge_gs}
\alias{merge_gs}


\title{
cluster/merge GatingSets based on the gating tree structures.
}
\description{
merge GatingSets based on the gating tree structures.
}
\usage{
merge_gs(x, ...)


}
\arguments{
  \item{x}{A \code{list} of \code{GatingSet}s .}
}
\details{
Group the individual GatingSets by their gating schemes.It is done by comparing the node list returned by \code{\link{getNodes}},which assumes
they follow the same population naming conventions.
Meanwhile the unused channels are automatically dropped to make sure the flow data has identical data structure within each group.
In order to further merge multiple GatingSet objects into one, use \code{\link{rbind2}}.     
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
    
	gs_list <- merge_gs(gs_list)
	gs_list
	
}
}
\keyword{ manip }
