\name{getGapdh3}
\alias{getGapdh3}

\title{ Gets the probesetid for the specified control probe }
\description{
  Gets the probesetid for the specified control probe
}
\usage{
  getGapdh3(name)
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
     probesetid <- getGapdh3(cleancdfname(cdfName((eset))))
  }
}

\keyword{ misc }
