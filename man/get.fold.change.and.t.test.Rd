\name{get.fold.change.and.t.test}
\alias{get.fold.change.and.t.test}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compute pairwise comparison statistics between two experimental groups }
\description{
  Generate fold changes (and possibly means) for a pair of experimental groups
}
\usage{
  get.fold.change.and.t.test(x,group,members,logged = T, return.exprs = "all",a.order=NULL,b.order=NULL,method=c("unlogged","logged","median"));
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an AffyBatch object }
  \item{group}{ column in pData(x) }
  \item{members}{ labels in group }
  \item{logged}{ is the AffyBatch data logged? }
  \item{return.exprs}{ Return the original expression levels for each chip. "all" returns for all arrays in the experiment, "samples", just for the pair being compared }
  \item{a.order}{ For a pairwise comparison the ordering of the first group of replicates}
  \item{b.order}{ For a pairwise comparison the ordering of the second group of replicates}
  \item{method}{ What method should be used to calculate the average for the fold-change - can be either "logged","unlogged","median" }
}
\details{
  Given an AffyBatch object, generate quick stats for pairwise comparisons between a pair of experimental groups. If a.order and b.order are specified then a paired sample t-test will be conducted between the groups, with the arrays in each group sorted according to the ordering specified. 
  The fold-changes are computed from the average values across replicates. By default this is done using the mean of the unlogged values. The parameter, method allows the mean of the logged values or the median to be used instead. T-tests are always computed with the logged data.
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
   eset.rma <- call.exprs(eset,"rma");
   results.a.v.none  <- get.fold.change.and.t.test(eset.rma,"drug",c("a","none"));
   results.a.v.none.pairwise  <- get.fold.change.and.t.test(eset.rma,"drug",c("a","none"),a.order=c(1,3,2),b.order=(1,2,3));
  }
}
\keyword{ misc }
