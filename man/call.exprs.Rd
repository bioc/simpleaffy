\name{call.exprs}
\alias{call.exprs}
\title{ normalise Affymetrix data }
\description{
  Generates expression summaries and normalises Affy data using either MAS5.0, GCRMA or RMA algorithms.
}
\usage{
call.exprs(x, algorithm = "rma", do.log = TRUE, sc = 100, method = NA)
}

\arguments{
  \item{x}{ an \code{AffyBatch} object }
  \item{algorithm}{ one of '"rma","rma-R","gcrma" or "mas5". "rma" makes use of a native C-library and  is faster than "rma-R" }
  \item{do.log}{ return logged data if true }
  \item{sc}{ if the mas5.0 algorithm is being used, sets the target intensity to scale to. Defaults to 100. }
  \item{method}{ The algorithm used to normalise the data. Has no effect for "rma", 
                 defaults to quantile normalisation for "rma" and no normalisation for "mas5" }
}
\details{

}
\value{
  An AffyBatch object containing expression summaries
}
\references{ http://bioinformatics.picr.man.ac.uk/ }
\author{ Crispin J Miller }

\seealso{ \code{\link{read.affy}}, \code{\link{expresso}}, \code{\link{rma}}}
\keyword{ misc }
\examples{
  \dontrun{
    eset.rma <- call.exprs(eset,"rma");
  }
}
