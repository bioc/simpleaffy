\name{getAlpha1}
\alias{getAlpha1}

\title{ Gets the probesetid for the specified control probe }
\description{
  Gets the probesetid for the specified control probe
}
\usage{
  getAlpha1(name)
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
     probesetid <- getAlpha1(cleancdfname(cdfName((eset))))
  }
}

\keyword{ misc }
