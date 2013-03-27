\name{getData}
\alias{getData}
\alias{getData-methods}
\alias{getData,GatingSetInternal,name-method}


\title{
Return the flowSet associated with a GatingSet by boolean expression
}
\description{
Return the flowSet associated with a GatingSet by boolean expression
}
\usage{
\S4method{getData}{GatingHierarchy}(obj,y,...)

}
\arguments{
  \item{obj}{
A \code{GatingSet} object . 
}
  \item{y}{
\code{name} boolean expression specifying the boolean combination of different cell populations
}
}
\details{
Returns a flowSet containing the events defined at by boolean expression \code{y}.
} 
\value{

A \code{ncdfFlowSet} or \code{flowSet}. 
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
\code{\link{getIndices}} \code{\link{getProp}} \code{\link{getPopStats}}
}

\examples{
\dontrun{
	#G is a GatingSet
	geData(G,3)
	res <- getData(gs[1],quote(`4+/TNFa+|4+/IL2+`))
	res[[1]]

	
}
}
\keyword{ manip }
