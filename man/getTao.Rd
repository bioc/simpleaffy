\name{getTao}
\alias{getTao}

\title{ Gets the Tao value the specified control probe }
\description{
  Gets the Tao value for the specified control probe
}
\usage{
  getTao(name)
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
     tao <- getTao(cleancdfname(cdfName((eset))))
  }
}

\keyword{ misc }
