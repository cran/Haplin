\name{plot.haplin}
\alias{plot.haplin}
%\alias{plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot a haplin object}
\description{
 Plot a haplin object and (optionally) produce picture files
}
\usage{
\method{plot}{haplin}(x, reference, separate.plots = F, filename, 
filetype = "png", use.dd, ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{}{Of the following arguments, only \code{x} is required.}
  \item{x}{A \code{haplin} object, i.e. the result of running \code{haplin}.}
  \item{reference}{Same as \code{reference} argument in \code{haplin}. Note that when plotting, you can only choose \code{"reciprocal"}, \code{"population"} or \code{"ref.cat"}. You cannot use a numeric value to change the reference category, to do that \code{haplin} must be run over again. (See the \code{reference} argument of \code{haplin}.)}
  \item{separate.plots}{Logical. If you estimate effects of both fetal and maternal genes you can decide whether or not to plot them in the same plot. The default is the same plot (TRUE), the alternative (FALSE) means in separate plots. If you choose separate plots you may have to set the graphics window to "recording" to make sure you can scroll back to the first plot.}
  \item{filename}{If you want a file containing the plot to be produced, give a character string for the filename.}
  \item{filetype}{The default filetype is \code{"png"}, alternatively you can choose \code{"jpeg"}.}
  \item{use.dd}{Numeric vector indicating which double dose estimates should be plotted. For instance, if set to c(1,3) only the first and third haplotypes will be drawn with double dose estimates. This is useful if some haplotypes are rare and you want to exclude the uncertain estimates from the plot.}
  \item{...}{Further arguments to be passed on to the plot function}
}
%\details{}
\references{ Gjessing HK and Lie RT. Case-parent triads: Estimating single- and double-dose effects of fetal and maternal disease gene haplotypes. Annals of Human Genetics (2006) 70, pp. 382-396.\cr\cr
Web Site: \url{http://folk.uib.no/gjessing/genetics/software/haplin/}}
\author{Hakon K. Gjessing\cr
Professor of Biostatistics\cr
Division of Epidemiology\cr
Norwegian Institute of Public Health\cr
\email{hakon.gjessing@uib.no}}
\note{Further information is found on the web page.
}
\seealso{\code{\link{haplin}}}

\examples{
\dontrun{

# Produce separate plots for child and mother, dump plots to files:
res <- haplin("data.dat", use.missing = T, maternal = T)
plot(res, separate.plots = T, filename = "Haplinres.png")

}}




