\name{getAlpha1}
\alias{getAlpha1}
\alias{getAlpha2}
\alias{setAlpha}

\title{ Gets the alpha1 or alpha2 values for the specified array }
\description{
  Gets (or sets) the alpha values for the specified control probe. Note
  that the name should match the results of \code{cleancdfname(cdfName((eset))} for the
  ExpressionSet of interest.
}
\usage{
  getAlpha1(name)
  getAlpha2(name)
  setAlpha(name,a1,a2)
}

\arguments{
  \item{name}{ A valid (clean) cdf name for an array }
  \item{a1}{ Alpha1 value}
  \item{a2}{ Alpha2 value}
}

\value{
  A probesetid
}

\author{ Crispin J Miller }

\examples{
  \dontrun{
     alpha1 <- getAlpha1(cleancdfname(cdfName((eset))))
  }
}

\keyword{ misc }
