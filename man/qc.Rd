\name{qc}
\alias{qc}
\alias{qc,AffyBatch-method}

\title{ Generate QC stats from an AffyBatch object }
\description{
  ...
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
