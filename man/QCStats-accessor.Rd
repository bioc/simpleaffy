\name{spikeInProbes}
\alias{spikeInProbes}
\alias{qcProbes}
\alias{avbg}
\alias{maxbg}
\alias{minbg}
\alias{percent.present}
\alias{sfs}
\alias{target}
\alias{ratios}
\title{ Accessor Methods for QCStats class }

\usage{
spikeInProbes(object)
qcProbes(object)
avbg(object)
maxbg(object)
minbg(object)
percent.present(object)
sfs(object)
target(object)
ratios(object)
}

\description{

Affymetrix recommend a series of QC metrics that should be used to
check that arrays have hybridised correctly and that sample quality is
acceptable. These are discussed in the document 'QC and Affymetrix
data' accompanying this package, and on the web at
http://bioinformatics.picr.man.ac.uk. They are described in detail in
the 'Expression Analysis Fundamentals' manual available from
Affymetrix.

Each of these functions queries a QCStats object to retrieve one of
these QC metrics. Refer to the accompanying documentation for more details
about the meaning of these values.

Briefly:


spikeInProbes:   Returns a table of intensities for the BioB, BioC, BioD and CreX probesets.
qcProbes:        Returns a table of intensities for GAPDH and BetaActin 3', 5' and Mid' probesets.
avbg:            Returns a vector of the average background levels for each array.
maxbg:           Returns a vector of the maximum background levels for each array.
minbg:           Returns a vector of the minimum background levels for each array.
percent.present: Returns a vector listing the percentage of probesets called present on each array.
sfs:             Returns a vector of scale factors for each array (as produced by the MAS 5.0 algorithm).
target:          The target intensity the arrays were scaled to.
ratios:		 A table of 3'/5' and 3'/M' ratios for GAPDH and beta-actin


}
\seealso{\code{\link{qc}}}  
\arguments{
  \item{object}{ the QCStats object you want to access }
}

\author{ Crispin J Miller }

\keyword{ misc }
