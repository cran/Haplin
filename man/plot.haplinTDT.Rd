\name{plot.haplinTDT}
\alias{plot.haplinTDT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot a haplinTDT object}
\description{
  Plots a haplinTDT object.
}
\usage{
  \method{plot}{haplinTDT}(x, separate.plots = FALSE, filename, filetype = "png", 
    ask = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  Of the following arguments, only 'x' is required.
  \item{x}{A \code{\link{haplinTDT}} object, i.e. the result of running
    \code{\link{haplinTDT}}.}
  \item{separate.plots}{Logical. The default is several plots on one page
    (FALSE), the alternative (TRUE) is one plot per page.}
  \item{filename}{If you want a file containing the plot to be produced, give a
    character string for the filename.}
  \item{filetype}{The default filetype is "png", alternatively you can
    choose "jpeg" or "pdf".}
  \item{ask}{Logical; if 'TRUE', the user is asked before each new page
    to be plotted.}
  \item{...}{Further arguments to be passed on to the plot function.}
}
%\references{}

\author{Hakon K. Gjessing (Professor of Biostatistics) and Oivind Skare\cr
  Division of Epidemiology\cr
  Norwegian Institute of Public Health\cr
  \email{hakon.gjessing@fhi.no}}


%\note{Further information is found on the web page.
%}
\seealso{\code{\link{haplinTDT}}, \code{\link{summary.haplinTDT}}}
\examples{

\dontrun{

## Standard run:
res <- haplinTDT("data.dat", nsim.perm=1000)
plot(res)
}
}