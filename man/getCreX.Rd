\name{getCreX}
\alias{getCreX}

\title{ Gets the probesetid for the specified control probe }
\description{
  Gets the probesetid for the specified control probe
}
\usage{
  getCreX(name)
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
     probesetid <- getCreX(cleancdfname(cdfName((eset))))
  }
}

\keyword{ misc }
