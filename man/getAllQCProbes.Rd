\name{getAllQCProbes}
\alias{getAllQCProbes}

\title{ Get the names of all the QC (not spike-in) probes on the array }
\description{
  What QC probes are there on the chip?
}
\usage{
  getAllQCProbes(name)
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
     qcprobes <- getAllQCProbes(cleancdfname(cdfName((eset))))
  }
}

\keyword{ misc }
