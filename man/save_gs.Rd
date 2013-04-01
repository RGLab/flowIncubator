\name{save_gs}
\alias{save_gs}
\alias{load_gs}

\title{
save/load a GatingSetInternal to/from disk.
}
\description{
Save/load a GatingSetInternal which is the gated flow data including gates and populations to/from the disk. 
The \code{GatingSetInternal} object The internal C data structure (gating tree),\code{ncdfFlowSet} object(if applicable) 
}
\usage{
save(G,path,overwrite = FALSE,...) 
load_gs(path)
}
\arguments{
  
  \item{G}{ A \code{GatingSetInternal}}
  
  \item{path}{A character scalar giving the path to save/load the GatingSet to/from.}

  \item{overwrite}{A logical scalar specifying whether to overwrite the existing folder.}
   
   \item{...}{other arguments:
   
   			 copy.cdf = TRUE: a logical scalar indicating whether to copy the cdf data file. Sometime it is more efficient to move the existing cdf file to the archived folder than to copy it. 
   							  however, by moving instead of copying the cdf file, the cdf is not accessible to the original GatingSet object anymore.
   							  So it is up to user to turn this "cp" off and manually do the "mv" the file to the destination folder.  
   				}    					
}

\value{
  
  \code{load_gs} returns a GatingSetInternal object
  
}

\author{
Mike Jiang


Maintainer: Mike Jiang <wjiang2@fhcrc.org>
}

\seealso{\code{\link{GatingSetInternal-class}}}

\examples{
\dontrun{
	#G is a GatingSet
	save_gs(G,path="tempFolder")
	G1<-load_gs(path="tempFolder")
}

}
\keyword{IO}

