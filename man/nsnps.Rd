% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsnps.R
\name{nsnps}
\alias{nsnps}
\title{Count the number of markers in the data}
\usage{
nsnps(data.in, design = "triad")
}
\arguments{
\item{data.in}{The data read in by \link{genDataRead}.}

\item{design}{The design used in the study - choose from:
\itemize{
  \item \emph{triad} (default) - data includes genotypes of mother, father and child;
  \item \emph{cc} - classical case-control;
  \item \emph{cc.triad} - hybrid design: triads with cases and controls
}}
}
\value{
How many markers (integer).
}
\description{
This is a help function to count the number of markers in an object read in with
 \link{genDataRead} (or loaded with \link{genDataLoad}).
}
