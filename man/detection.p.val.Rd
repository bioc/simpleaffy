\name{detection.p.val}
\alias{detection.p.val}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Filter pairwise comparison statistics between two experimental groups }
\description{
  Calculate detection p values and Present Marginal Absent calls, a la MAS5.0. This is an implementation based i nthe algorithm described in 
  Liu, Mei et al. (2002) 'Analysis of high density expression microarrays with signed-rank call algorimths', Bioinformatics 18(12) pp1593-1599.
}
\usage{
  detection.p.val(x, tao = 0.015,calls=T,alpha1=0.04,alpha2=0.06, ignore.saturated=T)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}              { An unnormalised AffyBatch object }
  \item{tao}            { Errrmmm... tao }
  \item{alpha1}         { Present-Marginal threshold }
  \item{alpha2}         { Marginal-Absent threshold }
  \item{calls}          { if true generate PMA calls }
  \item{ignore.saturated} { if true do the saturation correction described in the paper, with a saturation level of 46000 }
}
\value{
  A list containing "detection.p.values" and "calls".
}
\references{ http://bioinformatics.picr.man.ac.uk/ }
\author{ Crispin J Miller }

 ~Make other sections like Warning with \section{Warning }{....} ~ 

\examples{
 \dontrun{
   eset     <- read.affy();
   detection.p.val(eset);
  }
}
\keyword{ misc }
