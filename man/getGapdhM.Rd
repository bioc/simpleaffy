\name{getGapdhM}
\alias{getGapdhM}

\title{ Gets the probesetid for the specified control probe }
\description{
  Gets the probesetid for the specified control probe
}
\usage{
  getGapdhM(name)
}

\arguments{
   \item{name}{ A valid (clean) cdf name for an array }
}

\value{
  A probesetid
}

\author{ Crispin J Miller }

\examples{
  \dontrun{
     probesetid <- getGapdhM(cleancdfname(cdfName((eset))))
  }
}

\keyword{ misc }
