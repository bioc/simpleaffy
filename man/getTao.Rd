\name{getAlpha2}
\alias{getAlpha2}

\title{ Gets the probesetid for the specified control probe }
\description{
  Gets the probesetid for the specified control probe
}
\usage{
  getAlpha2(name)
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
     probesetid <- getAlpha2(cleancdfname(cdfName((eset))))
  }
}

\keyword{ misc }
