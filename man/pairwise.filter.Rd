\name{pairwise.filter}
\alias{pairwise.filter}

\title{ Filter pairwise comparison statistics between two experimental groups }

\description{
  Given the results of a pairwise.comparison, filter the resulting gene list on expression level, PMA calls (if available), fold change and t-test statistic
}

\usage{
  pairwise.filter(object,x,min.exp=log2(100),min.exp.no=0,min.present.no=3,fc=1.0,tt=0.001)
}


\arguments{
  \item{object}{ a 'PairComp' object }
  \item{x}{ an 'exprSet' object used to filter by intensity }
  \item{min.exp}{ Filter genes using a minimum expression cut off }
  \item{min.exp.no}{ A gene must have an expression intensity greater than 'min.exp' in at least this number of chips }
  \item{min.present.no}{ A gene must be called present on at least this number of chips }
  \item{fc}{ A gene must show a log2 fold change greater than this to be called significant }
  \item{tt}{ A gene must be changing with a p-score less than this to be called significant }
}
\value{
  A 'PairComp' object filtered to contain only the genes that pass the specified filter parameters.
}
\references{ http://bioinformatics.picr.man.ac.uk/ }
\author{ Crispin J Miller }



\examples{ 
  \dontrun{
    pc <- pairwise.comparison(eset.rma,"group",c("A","P"))
    pf <- pairwise.filter(pc,eset,tt=0.01);
  }
}

\keyword{ misc }
