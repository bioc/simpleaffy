\name{PairComp-class}
\docType{class}
\alias{PairComp-class}
\alias{[,PairComp-method}
\alias{[<-,PairComp-method}
\alias{calls,PairComp-method}
\alias{fc,PairComp-method}
\alias{group,PairComp-method}
\alias{means,PairComp-method}
\alias{members,PairComp-method}
\alias{pairwise.filter,PairComp-method}
\alias{tt,PairComp-method}
\alias{calls}
\alias{fc}
\alias{group}
\alias{means}
\alias{members}
\alias{tt}

\title{Class "PairComp" Represents the results of pairwise comparison between two experimental factors }
\description{  Holds fold-change, ttest p-score and detection p-value calls(if used) between a pair of experimental factors.}

\section{Slots}{
  \describe{
    \item{\code{means}:}{Object of class \code{"matrix"} Mean values for each of the experimental factors. }
    \item{\code{fc}:}{Object of class \code{"numeric"} Fold change between the means }
    \item{\code{tt}:}{Object of class \code{"numeric"} P-score between the factors }
    \item{\code{calls}:}{Object of class \code{"matrix"} Detection p-values for each probeset on each array }
    \item{\code{group}:}{Object of class \code{"character"} The name of the factor that was compared }
    \item{\code{members}:}{Object of class \code{"character"} A list containing the two levels compared between }
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "PairComp")}: get the values for the specified gene(s) }
    \item{[<-}{\code{signature(x = "PairComp")}: not supported. }
    \item{calls}{\code{signature(object = "PairComp")}: the detection.p.values }
    \item{fc}{\code{signature(object = "PairComp")}: the fold-changes }
    \item{group}{\code{signature(object = "PairComp")}: the name of the group that was compared }
    \item{means}{\code{signature(object = "PairComp")}: the means of the two experimental factors that were compared}
    \item{members}{\code{signature(object = "PairComp")}: the members of that group that were compared }
    \item{pairwise.filter}{\code{signature(object = "PairComp")}: Take a \code{PairComp} object and filter it to yield probesets that pass the specified criteria }
    \item{tt}{\code{signature(object = "PairComp")}: the results of a ttest between groups }
  }
}

\author{ Crispin Miller }

 ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{

}
\keyword{classes}