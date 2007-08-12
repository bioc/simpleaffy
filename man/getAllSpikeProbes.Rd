\name{getAllSpikeProbes}
\alias{getAllSpikeProbes}
\alias{setAllSpikeProbes}

\title{ Get the names of all the Spike (not QC) probes on the array }
\description{
  What Spike probes are there on the chip?
}
\usage{
  getAllSpikeProbes(name)
  setAllSpikeProbes(name,bioB,bioC,bioD,creX)

}

\arguments{
  \item{name}{ A valid (clean) cdf name for an array }
  \item{bioB}{bioB probeset name}
  \item{bioC}{bioC probeset name}
  \item{bioD}{bioD probeset name}
  \item{creX}{creX probeset name}
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
