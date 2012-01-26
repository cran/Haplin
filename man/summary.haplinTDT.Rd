\name{summary.haplinTDT}
\alias{summary.haplinTDT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summary of a haplinTDT object}
\description{
  Provides detailed information about a
  \code{haplinTDT} object.
}
\usage{
\method{summary}{haplinTDT}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A \code{haplinTDT} object, i.e. the result of running \code{haplinTDT}.}
  \item{...}{Further arguments to be passed on to the summary function.}
}
%\references{}

\author{Hakon K. Gjessing (Professor of Biostatistics) and Oivind Skare\cr
  Division of Epidemiology\cr
  Norwegian Institute of Public Health\cr
  \email{hakon.gjessing@fhi.no}}

%\note{}

\seealso{\code{\link{haplinTDT}}, \code{\link{plot.haplinTDT}}}
\examples{

\dontrun{

# Standard run with permutation test:
res <- haplinTDT("data.dat", nsim.perm=1000)
summary(res)
}
}