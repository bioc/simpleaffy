\name{getGapdh5}
\alias{getGapdh5}

\title{ Gets the probesetid for the specified control probe }
\description{
  Gets the probesetid for the specified control probe
}
\usage{
  getGapdh5(name)
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
     probesetid <- getGapdh5(cleancdfname(cdfName((eset))))
  }
}

\keyword{ misc }
