\name{read.affy}
\alias{read.affy}

\title{ Read a set of affy .CEL files and associated phenotypic data }
\description{
 Reads the file specified by covdesc, which defines phenotypic data for a
 set of .cel files. Reads the specified files into an AffyBatch object and
 then creates a pData object, defining the experimental factors 
 or covariants for those chips.
}
\usage{
read.affy(covdesc = "covdesc", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{covdesc}{ 
    A white space delimited file suitable for reading as a data frame. The first column (with no column name)  
     defines a list of .cel files to read, remaining columns (with names) represent experimental factors
     for each chip,
  }
  \item{\dots}{ extra functions to pass on to \code{ReadAffy} }
}
\details{ }
\value{
  An AffyBatch object
}
\references{ http://bioinformatics.picr.man.ac.uk/ }
\author{ Crispin J Miller }

\seealso{ \code{\link{ReadAffy}}, \code{\link{AffyBatch}} \code{\link{data.frame}}}
\keyword{ misc }
\examples { 
  \dontrun{
    eset <- read.affy();
    eset.rma <- call.exprs(eset,"rma");
  }
} 
