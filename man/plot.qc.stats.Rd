\name{plot.qc.stats}
\alias{plot.qc.stats}
\alias{plot,QCStats}
\alias{plot,QCStats-method}

\title{ Plots a QCStats object }
\description{
  Generates a visual summary of the various QC statistics recommended by Affymetrix in their 'Data Analysis Fundamentals' handbook.
}
\section{Usage}{
plot.qc.stats(x, fc.line.col = "black", sf.ok.region = "light blue", chip.label.col = "black", sf.thresh = 3, gdh.thresh = 3, ba.thresh = 3, present.thresh = 10, bg.thresh = 20, label = NULL, ...)
}
\arguments{
  \item{x}{ A \code{QCStats} object }
  \item{fc.line.col}{ The colour to mark fold change lines with }
  \item{sf.ok.region}{ The colour to mark the region in which scale factors lie within appropriate bounds }
  \item{chip.label.col}{ The colour to label the chips with }
  \item{sf.thresh}{ Scale factors must be within this fold-range }
  \item{gdh.thresh}{ Gapdh ratios must be within this range }
  \item{ba.thresh}{ beta actin must be within this range }
  \item{present.thresh}{ The percentage of genes called present must lie within this range }
  \item{bg.thresh}{ Array backgrounds must lie within this range }
  \item{label}{ What to call the chips }
  \item{\dots}{ Other parameters to pass through to \code{plot} }
}
\details{
  A lot of information is presented in this one figure. The dotted
circles correspond to fold-change - the outer one, labelled '+3',
corresponds to a 3-fold up-scaling, the
middle to no change, and the inner to a 3-fold
down-scaling. The blue circle represents the region in which all
arrays have scale factors within, by default, three-fold of each other. Its
position is found by calculating the mean scale factor for all chips
and positioning the center of the region so that the borders
are -1.5 fold up or down from the mean value.

Each array is plotted as a line from the 0-fold circle to the point
that corresponds to its scale factor. If the ends of all of the lines
are in the blue region, their scale-factors are compatible. The lines
are coloured blue if OK, red if not.

The figure also shows GAPDH and beta-actin 3'/5' ratios. These are
represented as a pair of points for each chip. Affy state that they
should be within 3. Any that fall outside the three-fold line are
coloured red; the rest are blue. 

Written around the edge of the figure are the number of genes called
present on each array and the average background. These will vary
according to the samples being processed, and Affy's QC suggests
simply that they should be similar. If any chips have significantly
different values this is flagged in red, otherwise the numbers are
displayed in blue. By default, 'significant' means that \%-present are
within 10\% of each other; background intensity, 20 units. These last
numbers are somewhat arbitrary and may need some tweaking to find
values that suit the samples you're dealing with, and the overall
nature of your Affymetrix setup.

Finally, if BioB is not present on a chip, this will be flagged by
printing 'BioB' in red next to the scale-factor for that chip.

In short, everything in the figure should be blue - red highlights a problem!
}
\author{ Crispin J Miller }

\seealso{ \code{\link{qc}} }
\examples{
  \dontrun{
    qc <- qc(eset,eset.mas)
    plot(qc)
  }
}


\keyword{ misc }
