\name{volcano}
\alias{volcano}

\title{ Draws a ttest v fold change plot from a PairComp object }
\description{
  Plota t-test pscore against fold change for a PairComp object. If a second object is supplied, identifies those points on the graph...
}
\usage{
  volcano(x,y=NULL,labels=colnames(means(x)))
}

\arguments{
  \item{x}{ a PairComp object to plot }
  \item{y}{ a PairComp object to highlight on the graph }
  \item{labels}{ used to generate the identifier on the fold-change axis}
}
\references{ http://bioinformatics.picr.man.ac.uk/ }
\author{ Crispin J Miller }

\seealso{ \code{\link{plot}}}


\examples{
  \dontrun{
    pc <- pairwise.comparison(eset.rma,"group")
    pcf <- pairwise.filter(eset.rma,eset,"group",tt=0.02,fc=1)

    volcano(pc,pcf)

  }
}

\keyword{ misc }
