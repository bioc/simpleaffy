\name{justMAS}
\alias{justMAS}


\title{ Generate Expression calls using a C implementation of the MAS 5.0 Algorithm}
\description{
  Implements the MAS5.0 background correction, expression summary and scaling functions as described in Affy's
  'Statistical Algorithms Description Document'
}
\usage{
  justMAS(unnormalised,tgt=100,scale=TRUE)
}


\arguments{
  \item{unnormalised}{ An unnormalised AffyBatch object }
  \item{tgt}{The target intensity to scale array to, if scaling.}
  \item{scale}{Scale the data to the specified target intensity.}
}

\details{ 
  Uses a C code implementation of the MAS5.0 algorithm (As described in Affymetrix's 'Statistical Algorithms Reference Guide' - see 
  \url{http://www.affymetrix.com/support/technical/technotes/statistical_reference_guide.pdf}).
}
\value{
  An AffyBatch object, with, in addition, scale-factors for each array stored in the object's \code{description@preprocessing@sfs} slot, and the target intensity the arrays were scaled to in  \code{description@preprocessing@tgt}
}

\url{ http://bioinformatics.picr.man.ac.uk/ }

\author{ Crispin J Miller }

\examples{ 
   \dontrun{
     eset.mas <- justMAS(eset.mas);    
   }
}

\seealso{\url{http://www.affymetrix.com/support/technical/technotes/statistical_reference_guide.pdf}}


\keyword{ misc }
