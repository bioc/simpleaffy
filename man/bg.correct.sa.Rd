\name{bg.correct.sa}
\alias{bg.correct.sa}


\title{ SimpleAffy implementation of mas5 background correction }
\description{
  Implements the MAS5.0 background correction functions as described in Affy's
  'Statistical Algorithms Description Document'. Should do the same thing as bg.correct.mas5.
}
\usage{
  bg.correct.sa(unnormalised,grid=c(4,4))
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{unnormalised}{ An unnormalised AffyBatch object }
  \item{grid}{The dimensions of the grid to divide the chip into for background correction.}
}
\value{
  An AffyBatch object
}
\references{ http://bioinformatics.picr.man.ac.uk/ }
\author{ Crispin J Miller }

 ~Make other sections like Warning with \section{Warning }{....} ~ 

\examples{
 \dontrun{
   eset     <- read.affy();
   justMas(eset);
  }
}
\keyword{ misc }
