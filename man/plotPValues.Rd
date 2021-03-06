% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotPValues.R
\name{plotPValues}
\alias{plotPValues}
\title{Plotting p-values for relative risks}
\usage{
plotPValues(
  object,
  windows,
  which.p.val = "overall",
  plot.signif.only = FALSE,
  signif.thresh = 0.05,
  title,
  filename
)
}
\arguments{
\item{object}{The \code{\link{haplinSlide}} results: list of \code{\link{haptable}}
objects.}

\item{windows}{Numerical vector; if given, the plot will be restricted to only those.}

\item{which.p.val}{Character string specifying which p-values to choose for plotting:
"overall" (default), "child", "child.double", "maternal", "maternal.double",
"paternal". The last three options can be chosen only if \code{\link{haplinSlide}}
 was run with \code{maternal = TRUE} or \code{poo = TRUE}.}

\item{plot.signif.only}{Logical: whether to filter out the "non-significant" markers
from the plot. Default: FALSE, i.e., plot everything.}

\item{signif.thresh}{The threshold defining the significant p-values: if
\code{plot.signif.only == TRUE}, then only the markers with relative risk p-values
lower than the threshold will be kept for plotting. Default: 0.05.}

\item{title}{Optional character string for the title of the figure.}

\item{filename}{If the plot should be saved to the disk, give the name of the output
file including the file extension.}
}
\value{
Invisibly returns the table with only the plotted p-values.
}
\description{
This function plots p-values for the relative risks calculated by \code{\link{haplinSlide}}.
}
\details{
The output of \code{\link{haplinSlide}} can be very lengthy and not suitable for
  an overall plot of all the relative risks (RR) on one figure. Therefore, it's advised
  to first plot only the p-values for each window (user can choose which p-values to
  plot - see parameter \code{which.p.val}), and only then plot the RRs for specific
  windows, for which the p-values are significant.
}
