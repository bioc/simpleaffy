\name{qc.affy}
\alias{qc.affy}

\title{ Generate Affymetrix Style QC Statistics }
\description{
  Generates Affymetrix style QC statistics as defined in Affymetrix' 'Data Analysis Fundamentals' handbook }
}
\usage{
  qc.affy(unnormalised, normalised = NULL, logged = TRUE, tau =
                 NULL, alpha1 = NULL, alpha2 = NULL, bioB = NULL, bioC
                 = NULL, bioD = NULL, creX = NULL, gapdh3 = NULL,
                 gapdhM = NULL, gapdh5 = NULL, actin3 = NULL, actinM =
                 NULL, actin5 = NULL)
}

\arguments{
  \item{unnormalised}{ An unnormalised raw \code{AffyBatch} object to call qc stats on }
  \item{normalised}{ The same one, processed using \code{\link{justMAS}} (contains scale factors etc.). If not supplied, then the object gets calculated internally. }
  \item{logged}{ True if the data is logged }
  \item{tau}{ Specify the tau value rather than the default for the arrays}
  \item{alpha1}{Specify the alpha1 value rather than the default for the arrays}
  \item{alpha2}{Specify the alpha2 value rather than the default for the arrays}
  \item{bioB}{Specify the bioB value rather than the default for the arrays}  
  \item{bioC}{Specify the bioC value rather than the default for the arrays}
  \item{bioD}{Specify the bioD value rather than the default for the arrays}  
  \item{creX}{Specify the crex value rather than the default for the arrays}  
  \item{gapdh3}{Specify the gapdh3 value rather than the default for the arrays}  
  \item{gapdhM}{Specify the gapdhM value rather than the default for the arrays}  
  \item{gapdh5}{Specify the gapdh5 value rather than the default for the arrays}  
  \item{actin3}{Specify the actin3 value rather than the default for the arrays}  
  \item{actinM}{Specify the actinM value rather than the default for the arrays}  
  \item{actin5}{Specify the actin5 value rather than the default for the arrays}
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
