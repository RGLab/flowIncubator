\name{merge}
\alias{merge}
\alias{split}
\alias{merge-methods}
\alias{split-methods}
\alias{split,GatingSetList-methods}
\alias{merge,GatingSetList-method}


\title{
split/merge GatingSets based on the gating tree structures.
}
\description{
split/merge GatingSets based on the gating tree structures.
}
\usage{
\S4method{merge}{GatingSetList}(x,path=tempdir(),...)
\S4method{split}{GatingSetList}(x,f="missing",...)

}
\arguments{
  \item{x}{
A \code{GatingSetList} object . 
}
\item{f}{not used}
  \item{path}{
\code{character} path to store the temporary cdf files 
}
}
\details{
split the individual GatingSets by their gating schemes.It is accomplished by comparing the node list returned by \code{\link{getNodes}},which assumes
they follow the same population naming conventions.
Merge does the splitting first and \code{\link[flowWorkspace::rbind2]{rbind2}} subsequently.Meanwhile the unused channels are automatically dropped to make sure the 
flow data has identical data structure.  
 }
\value{
\code{split} return a \code{list}.
\code{merge} return A \code{\link{GatingSetList}} that contains multiple GatingSets each of which share the same gating and data structure.
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
    #only splitting
	split(GatingSetList(gs_list))

	#splitting and merging
	gs_list_merged <- merge(GatingSetList(gs_list))
	#print the merged result
	gs_list_merged


	
}
}
\keyword{ manip }
