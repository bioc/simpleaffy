\name{trad.scatter.plot}
\alias{trad.scatter.plot}

\title{ Does a traditional scatter plot of expression data }
\description{
  Plots expression data as a scatter plot with optional fold-change lines
}
\usage{
trad.scatter.plot(x, y, add = F, fc.lines = c(1, 2, 3, 4), draw.fc.lines = T, draw.fc.line.labels = T, fc.line.col = "lightgrey", pch = 20, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ x coords }
  \item{y}{ y coords }
  \item{add}{ add this data to an existing graph }
  \item{fc.lines}{ Vector of intervals at which to draw fold-change lines }
  \item{draw.fc.lines}{ Draw fold change lines? }
  \item{draw.fc.line.labels}{ Label the fold change lines with the fold changes they represent? }
  \item{fc.line.col}{ The colour to draw fold change lines }
  \item{pch}{ Plotting character to use for the scatter data (see \code{plot} for more details) }
  \item{\dots}{ Additional parameters to pass through to the underlying \code{plot} function }
}
\references{ http://bioinformatics.picr.man.ac.uk/ }
\author{ Crispin J Miller }

\seealso{ \code{\link{plot}},\code{\link{simple.mva.plot}} }
\examples{
  \dontrun{
   eset <- read.affy();
   rma  <- norm(eset,"rma" )    
   results <- pairwise.comparison(rma,"drug",c("a","b"))

  # plot four scatter plots with differing p-value cutoffs
  # on each one show the effect of filtering by expression level by colouring
  # spots black,red,orange or yellow according to their intensity

  dev.off();
  cols=c("red","orange","yellow");
  par(mfrow=c(2,2));
  for(tt.cutoff in c(0.05,0.01,0.005,0.001)) {
    trad.scatter.plot(results[,"a"],results[,"b"]);
    legend(4,14,legend=paste("I >",c("0","10","100","200")),pch=20,col=c("black",cols));
    text(12,4,paste("p <",tt.cutoff));
    idx = 1;
    for(intensity.cutoff in log2(c(10,100,200))) {
      sig     <- pairwise.filter(results,min.exp=intensity.cutoff,fc=fc.cutoff,tt=tt.cutoff);
      trad.scatter.plot(sig[,"a"],sig.mc.rma[,"b"],add=T,col=cols[idx]);
      idx <- idx + 1;
    }
  }
  }  
}  
\keyword{ misc }
