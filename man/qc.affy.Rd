\name{qc.affy}
\alias{qc.affy}

\title{ Generate Affymetrix Style QC Statistics }
\description{
  Generates Affymetrix style QC statistics as defined in Affymetrix' 'Data Analysis Fundamentals' handbook }
}
\usage{
  qc.affy(unnormalised, normalised = NULL, logged=TRUE)
}

\arguments{
  \item{unnormalised}{ An unnormalised raw \code{AffyBatch} object to call qc stats on }
  \item{normalised}{ The same one, processed using \code{\link{justMAS}} (contains scale factors etc.). If not supplied, then the object gets calculated internally. }
  \item{logged}{ True if the data is logged }
}

\value{
  A QCStats object describing the supplied \code{\link{AffyBatch}}
}

\author{ Crispin J Miller }

\examples{
  \dontrun{
     qc <- qc(eset,eset.mas)
  }
}

\keyword{ misc }
