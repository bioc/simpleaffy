\name{getActin3}
\alias{getActin3}

\title{ Gets the probesetid for the specified control probe }
\description{
  Gets the probesetid for the specified control probe
}
\usage{
  getActin3(name)
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
     probesetid <- getActin3(cleancdfname(cdfName((eset))))
  }
}

\keyword{ misc }
