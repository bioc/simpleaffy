\name{QCStats-class}
\docType{class}
\alias{QCStats-class}

\alias{avbg,QCStats-method}
\alias{maxbg,QCStats-method}
\alias{minbg,QCStats-method}
\alias{percent.present,QCStats-method}
\alias{sfs,QCStats-method}
\alias{target,QCStats-method}
\alias{spikeInProbes,QCStats-method}
\alias{qcProbes,QCStats-method}
\alias{ratios,QCStats-method}

\alias{avbg-methods}
\alias{maxbg-methods}
\alias{minbg-methods}
\alias{percent.present-methods}
\alias{sfs-methods}
\alias{target-methods}
\alias{spikeInProbes-methods}
\alias{qcProbes-methods}
\alias{ratios-methods}

\title{Class "QCStats" }
\description{ Holds Quality Control data for a set of Affymetrix arrays }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{qc(AffyBatch)}.
}
\section{Slots}{
  \describe{
    \item{\code{scale.factors}:}{Object of class \code{"numeric"} Scale factors used to scale the chips to the specified target intensity }
    \item{\code{target}:}{Object of class \code{"numeric"} The target intensity to which the chips were scaled }
    \item{\code{percent.present}:}{Object of class \code{"numeric"} Number of genes called present }
    \item{\code{average.background}:}{Object of class \code{"numeric"} The average background for the arrays }
    \item{\code{minimum.background}:}{Object of class \code{"numeric"} The minimum background for the arrays }
    \item{\code{maximum.background}:}{Object of class \code{"numeric"} The maximum background for the arrays }
    \item{\code{spikes}:}{Object of class \code{"list"} spiked in probes (e.g. biob, bioc...) }
    \item{\code{qc.probes}:}{Object of class \code{"list"} qc probes (e.g. gapdh 3,5,M,...) }
  }
}
\section{Methods}{
  \describe{
    \item{avbg}{\code{signature(object = "QCStats")}:  average background}
    \item{maxbg}{\code{signature(object = "QCStats")}: maximum background }
    \item{minbg}{\code{signature(object = "QCStats")}: minimum background }
    \item{spikeInProbes}{\code{signature(object = "QCStats")}: the spike-in QC probes }
    \item{qcProbes}{\code{signature(object = "QCStats")}: the gapdh and actin QC probes }
    \item{percent.present}{\code{signature(object = "QCStats")}: no probesets called present }
    \item{plot}{\code{signature(x = "QCStats")}: Plot a QCStats object }
    \item{sfs}{\code{signature(object = "QCStats")}: scale factors }
    \item{target}{\code{signature(object = "QCStats")}: target scaling }
    \item{ratios}{\code{signature(object = "QCStats")}: 5'3' and 5'M ratios for QC Probes }
  }
}

\author{ Crispin J Miller }

\seealso{

}

\keyword{classes}
