\name{get.array.subset}
\alias{get.array.subset}

\title{ Get a subset of chips from an AffyBatch object according to their phenotypic data }
\description{
  Each AffyBatch object has a pData object, which is a data.frame, associated with it. Each column in the 
  data.frame (such as "group"), defines an experimental condition (e.g. "treatment" or "control"). This function
  extracts chips from an AffyBatch object that are in one or more category.
}
\usage{
get.array.subset(x, group, members)
}

\arguments{
  \item{x}{ An \code{AffyBatch} object }
  \item{group}{ A column name in the pData object }
  \item{members}{ a vector of labels within the specified column }
}
\details{
  given an AffyBatch object eset, with a column called "drug" containing labels "a", "b" and "none" in its pData, 
  get.array.subset(eset,"drug",c("a","none")) returns an AffyBatch object just containing the chips for "a" and "none".
}
\value{
  An AffyBatch object
}
\references{ http://bioinformatics.picr.man.ac.uk/ }
\author{ Crispin J Miller }
\keyword{ misc }
\examples{
  \dontrun{
    a.v.b <- get.array.subset(eset,"drug",c("a","b"));
  } 
}
