\name{QCStats-class}
\docType{class}
\alias{QCStats-class}
\alias{actin35,QCStats-method}
\alias{actin3,QCStats-method}
\alias{actin3M,QCStats-method}
\alias{actin5,QCStats-method}
\alias{actinM,QCStats-method}
\alias{avbg,QCStats-method}
\alias{bioB,QCStats-method}
\alias{bioC,QCStats-method}
\alias{bioD,QCStats-method}
\alias{creX,QCStats-method}
\alias{gapdh35,QCStats-method}
\alias{gapdh3,QCStats-method}
\alias{gapdh3M,QCStats-method}
\alias{gapdh5,QCStats-method}
\alias{gapdhM,QCStats-method}
\alias{maxbg,QCStats-method}
\alias{minbg,QCStats-method}
\alias{percent.present,QCStats-method}
\alias{sfs,QCStats-method}
\alias{target,QCStats-method}
\alias{actin35-methods}
\alias{actin3-methods}
\alias{actin3M-methods}
\alias{actin5-methods}
\alias{actinM-methods}
\alias{avbg-methods}
\alias{bioB-methods}
\alias{bioC-methods}
\alias{bioD-methods}
\alias{creX-methods}
\alias{gapdh35-methods}
\alias{gapdh3-methods}
\alias{gapdh3M-methods}
\alias{gapdh5-methods}
\alias{gapdhM-methods}
\alias{maxbg-methods}
\alias{minbg-methods}
\alias{percent.present-methods}
\alias{sfs-methods}
\alias{target-methods}

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
    \item{\code{bioB}:}{Object of class \code{"character"} bioB present call }
    \item{\code{bioC}:}{Object of class \code{"character"} bioC present call }
    \item{\code{bioD}:}{Object of class \code{"character"} bioD present call }
    \item{\code{creX}:}{Object of class \code{"character"} creX present call }
    \item{\code{gapdh3}:}{Object of class \code{"numeric"} gapdh 3' intensity }
    \item{\code{gapdhM}:}{Object of class \code{"numeric"} gapdh M  intensity }
    \item{\code{gapdh5}:}{Object of class \code{"numeric"} gapdh 5' intensity }
    \item{\code{actin3}:}{Object of class \code{"numeric"} beta actin 3' intensity }
    \item{\code{actinM}:}{Object of class \code{"numeric"} beta actin M  intensity }
    \item{\code{actin5}:}{Object of class \code{"numeric"} beta actin 5' intensity }
  }
}
\section{Methods}{
  \describe{
    \item{actin35}{\code{signature(object = "QCStats")}: beta actin 3' 5' ratio }
    \item{actin3}{\code{signature(object = "QCStats")}: beta actin 3' intensity }
    \item{actin3M}{\code{signature(object = "QCStats")}: beta actin 3' M ratio }
    \item{actin5}{\code{signature(object = "QCStats")}: beta actin 5' intensity }
    \item{actinM}{\code{signature(object = "QCStats")}: beta actin M intensity }
    \item{avbg}{\code{signature(object = "QCStats")}:  average background}
    \item{bioB}{\code{signature(object = "QCStats")}: bioB call }
    \item{bioC}{\code{signature(object = "QCStats")}: bioC call }
    \item{bioD}{\code{signature(object = "QCStats")}: bioD call }
    \item{creX}{\code{signature(object = "QCStats")}: creX call }
    \item{gapdh35}{\code{signature(object = "QCStats")}: gapdh 3' 5' ratio }
    \item{gapdh3}{\code{signature(object = "QCStats")}: gapdh 3' intensity }
    \item{gapdh3M}{\code{signature(object = "QCStats")}: gapdh 3' M ratio }
    \item{gapdh5}{\code{signature(object = "QCStats")}: gapdh 5' inensity }
    \item{gapdhM}{\code{signature(object = "QCStats")}: gapdh M  intensity }
    \item{maxbg}{\code{signature(object = "QCStats")}: maximum background }
    \item{minbg}{\code{signature(object = "QCStats")}: minimum background }
    \item{percent.present}{\code{signature(object = "QCStats")}: no probesets called present }
    \item{plot}{\code{signature(x = "QCStats")}: Plot a QCStats object }
    \item{sfs}{\code{signature(object = "QCStats")}: scale factors }
    \item{target}{\code{signature(object = "QCStats")}: target scaling }
  }
}

\author{ Crispin J Miller }

\seealso{

}

\keyword{classes}
