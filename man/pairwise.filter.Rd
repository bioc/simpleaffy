\name{pairwise.filter}
\alias{pairwise.filter}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Filter pairwise comparison statistics between two experimental groups }
\description{
  Given the results of a pairwise.comparison, filter the resulting gene list on expression level, PMA calls (if available), fold change and t-test statistic
}
\usage{
  pairwise.filter(pw,min.exp=log2(100),min.exp.no=0,min.present.no=3,fc=1.0,tt=0.001)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pw}{ a list resulting from the pairwise.comparison function }
  \item{min.exp}{ Filter genes using a minimum expression cut off }
  \item{min.exp.no}{ A gene must have an expression intensity greater than 'min.exp' in at least this number of chips }
  \item{min.present.no}{ A gene must be called present on at least this number of chips }
  \item{fc}{ A gene must show a log2 fold change greater than this to be called significant }
  \item{tt}{ A gene must be changing with a p-score less than this to be called significant }
}
\value{
  A list containing fold changes means t-tests detection p-values and PMA calls for each significant gene.
}
\references{ http://bioinformatics.picr.man.ac.uk/ }
\author{ Crispin J Miller }

 ~Make other sections like Warning with \section{Warning }{....} ~ 

\examples{
 \dontrun{
   eset     <- read.affy();
   eset.rma <- norm(eset,"rma");
   results.a.v.none  <- pairwise.comparison(eset.rma,"drug",c("a","none"));
   # find genes expressed with an intensity greater than 60 in at least 3 chips
   # with a fold change of more than 2 and a p-score less than 0.001
   significant.a.v.none  <- pairwise.filter(results.a.v.none,minexp=log2(60),min.no=3); 
   trad.scatter.plot(results.a.v.none$means[,"a"],results.a.v.none$means[,"none"]);
   trad.scatter.plot(results.a.v.none$means[names(significant.a.v.none$fc),"a"],results.a.v.none$means[names(significant.a.v.none$fc),"none"],add=T,col="red");
   anno <- get.annotation(names(significant.a.v.none$fc),"hgu133a");
   write.table(file="spreadsheet.xls",spreadsheet);
  }
}
\keyword{ misc }
