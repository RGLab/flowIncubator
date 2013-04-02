\name{merge}
\alias{merge}
\alias{merge-methods}
\alias{merge,GatingSetList-method}


\title{
cluster/merge GatingSets based on the gating tree structures.
}
\description{
split/merge GatingSets based on the gating tree structures.
}
\usage{
\S4method{merge}{GatingSetList}(x,path=tempdir(), test.run = TRUE, ...)


}
\arguments{
  \item{x}{
A \code{GatingSetList} object . 
}
  \item{path}{
				\code{character} path to store the temporary cdf files 
				}
 \item{test.run}{
			\code{logical} scalar indicating whether to perform the actual merging. 
			}
}
\details{
Group the individual GatingSets by their gating schemes.It is accomplished by comparing the node list returned by \code{\link{getNodes}},which assumes
they follow the same population naming conventions.
Merge does the splitting first and \code{\link[flowWorkspace::rbind2]{rbind2}} subsequently.Meanwhile the unused channels are automatically dropped to make sure the 
flow data has identical data structure.  
 }
\value{
A \code{list} that simply divides the original gatingSets into groups without actual merging the data when rest.run is TRUE.
When rest.run is FALSE, The merging is done and returns a \code{\link{GatingSetList}} that contains multiple GatingSets each of which share the same gating and data structure.
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
    #visual check the groups before the actual merging
	merge(GatingSetList(gs_list), test.run = TRUE)

	#do the merging
	gs_list_merged <- merge(GatingSetList(gs_list) , test.run = FALSE)
	#print the merged result
	gs_list_merged


	
}
}
\keyword{ manip }
