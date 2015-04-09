\name{getIndices}
\alias{getIndices}
\title{routine to return the indices by specify boolean combination of reference nodes:}
\arguments{
  \item{y}{a quoted expression.}
}
\description{
  It adds the boolean gates and does the gating on the fly,
  and return the indices associated with that bool gate,
  and remove the bool gate the typical use case would be
  extracting any-cytokine-expressed cells
}
\examples{
getIndices(gs,quote(`4+/TNFa+|4+/IL2+`)) (it may be faster than R version)
}

