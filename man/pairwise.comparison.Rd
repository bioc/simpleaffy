\name{pairwise.comparison}
\alias{pairwise.comparison}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compute pairwise comparison statistics between two experimental groups }
\description{
  Generate fold changes, t-tests and means for a pair of experimental groups
}
\usage{
  pairwise.comparison(x,group,members,spots)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an AffyBatch object }
  \item{group}{ column in pData(x) }
  \item{members}{ labels in group }
  \item{spots}{ unnormalised AffyBatch data for this experiment - if included, results in PMA calls and detection p-values being generated }
}
\details{
  Given an AffyBatch object, generate quick stats for pairwise comparisons between a pair of experimental groups.
  
}
\value{
  A list containing "fc" fold-changes, "means", means, "tt", t-test
p-values, and if 'spots' is present, "detection.p.values" detection p
values a la MAS5.0 and "calls" Present Marginal Absent calls.  

}
\references{ http://bioinformatics.picr.man.ac.uk/ } \author{ Crispin J Miller }

 ~Make other sections like Warning with \section{Warning }{....} ~ 

\examples{
 \dontrun{
   eset     <- read.affy();
   eset.rma <- norm(eset,"rma");
   results.a.v.none  <- pairwise.comparison(eset.rma,"drug",c("a","none"));
  }
}
\keyword{ misc }
