\name{get.fold.change.and.t.test}
\alias{get.fold.change.and.t.test}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compute pairwise comparison statistics between two experimental groups }
\description{
  Generate fold changes (and possibly means) for a pair of experimental groups
}
\usage{
  get.fold.change.and.t.test(x,group,members,logged = T, return.exprs = "all");
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an AffyBatch object }
  \item{group}{ column in pData(x) }
  \item{members}{ labels in group }
  \item{logged}{ is the AffyBatch data logged? }
  \item{return.exprs}{ Return the original expression levels for each chip. "all" returns for all arrays in the experiment, "samples", just for the pair being compared }
}
\details{
  Given an AffyBatch object, generate quick stats for pairwise comparisons between a pair of experimental groups.
  
}
\value{
  A list containing: "fc" fold changes,"tt", t-tests and "means", means for each gene.
}
\references{ http://bioinformatics.picr.man.ac.uk/ }
\author{ Crispin J Miller }

 ~Make other sections like Warning with \section{Warning }{....} ~ 

\examples{
 \dontrun{
   eset     <- read.affy();
   eset.rma <- norm(eset,"rma");
   results.a.v.none  <- get.fold.change.and.t.test(eset.rma,"drug",c("a","none"));
  }
}
\keyword{ misc }
