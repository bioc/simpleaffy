\name{read.affy}
\alias{read.affy}

\title{ Read a Set of .CEL Files and Phenotypic Data }
\description{
 Reads the specified file, which defines phenotypic data for a
 set of .CEL files. Reads the specified files into an \code{\link{AffyBatch}} object and
 then creates a \code{\link{phenoData}} object, defining the experimental factors 
 for those chips.
}
\usage{
read.affy(covdesc = "covdesc",path=".", ...)
}

\arguments{
  \item{covdesc}{ 
    A white space delimited file suitable for reading as a \code{\link{data.frame}}. The first column (with no column name)  
     contains the names(or paths to) the .CEL files to read. Remaining columns (with names) represent experimental factors
     for each chip. these become elements of the \code{\link{phenoData}} object.
  }
  \item{\dots}{ extra functions to pass on to \code{ReadAffy} }
  \item{path}{ The path to prefix the filenames with before calling \code{ReadAffy} }
}
\details{ }
\value{
  An AffyBatch object
}
\url{ http://bioinformatics.picr.man.ac.uk/ }
\author{ Crispin J Miller }

\seealso{ \code{\link{ReadAffy}}, \code{\link{AffyBatch}} \code{\link{data.frame}} \code{\link{phenoData}}}
\keyword{ misc }
\examples { 
  \dontrun{
    eset <- read.affy(); # read a set of CEL files
  }
  eset.rma <- call.exprs(eset,"rma");   
}
