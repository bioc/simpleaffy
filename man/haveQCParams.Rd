\name{haveQCParams}
\alias{haveQCParams}

\title{ True if data exists in the QC environment describing the specified chip }
\description{
  Looksin .qcEnv to see if the specified array is known about
}
\usage{
  haveQCParams(name)
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
     haveQCParams(cleancdfname(cdfName((eset))))
  }
}

\keyword{ misc }
