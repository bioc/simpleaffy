\name{qc.probes}
\docType{data}
\alias{qc.probes}
\title{ The name of the QC probes }
\description{
Contains the names of QC probes used on different chips
}
\details{
Affymetrix arrays contain a set of QC probes designed to check for
full length transcripts. Most arrays contain probes targeted at the 3'
5' and mid regions of the GAPDH and beta-actin transcripts. If all is
well, the ratios between signal at either end of the transcripts
should not be 'too large'. 

This object contains a table of the probeset ids used on each array.

For more details, see the document 'QC and
Affymetrix data' distributed with this package and available on the
web at http://bioinformatics.picr.man.ac.uk.


}
\seealso{\code{\link{getAllQCProbes}}, \code{\link{qcProbes}} }
\format{A CSV file, with a column for each spike.}
\keyword{datasets}
