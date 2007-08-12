\name{getAllQCProbes}
\alias{getAllQCProbes}
\alias{setAllQCProbes}

\title{ Get (or set) the names of all the QC (not spike-in) probes on the array }
\description{
  What QC probes are there on the chip?
  Note
  that the name should match the results of \code{cleancdfname(cdfName((eset))} for the
  ExpressionSet of interest.
}
\usage{
  getAllQCProbes(name)
  setAllQCProbes(name,a3,aM,a5,g3,gM,g5,rpol3=NA,rpolM=NA,rpol5=NA,tatabp3=NA,tatabpM=NA,tatabp5=NA)
}

\arguments{
  \item{name}{ A valid (clean) cdf name for an array }
  \item{a3}{actin 3' probeset name}
  \item{aM}{actin mid probeset name}
  \item{a5}{actin 5' probeset name}
  \item{g3}{gapdh 3' probeset name}
  \item{gM}{gapdh mid probeset name}
  \item{g5}{gapdh 5' probeset name}
  \item{rpol3}{RNA polymerase 3' probeset name}
  \item{rpolM}{RNA polymerase mid probeset name}
  \item{rpol5}{RNA polymerase 5' probeset name}
  \item{tatabp3}{tata binding protein 3' probeset name}
  \item{tatabpM}{tata binding protein mid probeset name}
  \item{tatabp5}{tata binding protein 5' probeset name}
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
