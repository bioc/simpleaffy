\name{qc}
\alias{qc}
\alias{qc,AffyBatch-method}

\title{ Generate QC stats from an AffyBatch object }

\description{
Generate QC metrix for Affymetrix data.
}

\details{
Affymetrix recommend a series of QC metrics that should be used to
check that arrays have hybridised correctly and that sample quality is
acceptable. These are discussed in the document 'QC and Affymetrix
data' accompanying this package, and on the web at
http://bioinformatics.picr.man.ac.uk. They are described in detail in
the 'Expression Analysis Fundamentals' manual available from
Affymetrix.

Before using this function you are strongly encouraged to read the 'QC
and Affymetrix data' doument, which contains detailed examples.

This function takes an \code{AffyBatch} object and generates a \code{QCStats}
object containing a set of QC metrics. See \code{qc.affy} for more details.
}
\usage{
qc(unnormalised, ...)
}
\arguments{
  \item{unnormalised}{ An AffyBatch object with nowt done to it }
  \item{\dots}{ Any other parameters }
}

\author{ Crispin J Miller }

\seealso{ \code{\link{qc.affy}} }
\examples{
  \dontrun{
    qc <- qc(eset,eset.mas)
    plot(qc)
  }
}


\keyword{ misc }
