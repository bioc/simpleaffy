\name{getAllSpikeProbes}
\alias{getAllSpikeProbes}

\title{ Get the names of all the Spike (not QC) probes on the array }
\description{
  What Spike probes are there on the chip?
}
\usage{
  getAllSpikeProbes(name)
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
     qcprobes <- getAllSpikeProbes(cleancdfname(cdfName((eset))))
  }
}

\keyword{ misc }
