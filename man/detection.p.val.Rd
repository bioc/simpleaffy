\name{detection.p.val}
\alias{detection.p.val}

\title{ Calculate Detection P Values  }
\description{
  Calculate detection p values and Present Marginal Absent calls, a la MAS5.0. This is an implementation based on the algorithm described in 
  Liu, Mei et al. (2002) 'Analysis of high density expression microarrays with signed-rank call algorimths', Bioinformatics 18(12) pp1593-1599.
}
\usage{
  detection.p.val(x, tau = 0.015,calls=TRUE,alpha1=0.04,alpha2=0.06, ignore.saturated=TRUE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ An unnormalised AffyBatch object }
  \item{tau}{ Errrmmm... tau }
  \item{alpha1}{ Present-Marginal threshold }
  \item{alpha2}{ Marginal-Absent threshold }
  \item{calls}{ if true generate PMA calls }
  \item{ignore.saturated}{ if true do the saturation correction described in the paper, with a saturation level of 46000 }
}
\value{
  A list:
  \item{pval}{A matrix of detection p values}
  \item{calls}{PMA calls}
}
\references{ http://bioinformatics.picr.man.ac.uk/ }
\author{ Crispin J Miller }

 ~Make other sections like Warning with \section{Warning }{....} ~ 

\note{
  Alpha1 and Alpha2 are parameters that change according to the chip type you are using. see \code{\link{getQCParams}} for more details.
}
\seealso{
  \code{\link{getQCParams}}
}
\examples{
   \dontrun{
     dpv <- detection.p.val(eset);
   }
  }
}
\keyword{ misc }
