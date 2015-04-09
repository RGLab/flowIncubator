\name{getData}
\alias{getData}
\title{Return the flowSet associated with a GatingSet by boolean expression}
\arguments{
  \item{obj}{A \code{GatingSet} object .}

  \item{y}{\code{name} boolean expression specifying the
  boolean combination of different cell populations}
}
\value{
  A \code{list} of \code{numerci matrices}
}
\description{
  Returns a flowSet containing the events defined at by
  boolean expression \code{y}.
}
\examples{
\dontrun{
	#G is a GatingSet
	geData(G,3)
	res <- getData(gs[1],quote(`4+/TNFa+|4+/IL2+`))
	res[[1]]
}
}
\author{
  Mike Jiang \email{wjiang2@fhcrc.org}
}
\seealso{
  \code{\link{getIndices}} \code{\link{getProp}}
  \code{\link{getPopStats}}
}

