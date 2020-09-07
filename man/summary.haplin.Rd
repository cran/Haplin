\name{summary.haplin}
\alias{summary.haplin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summary of a haplin object}
\description{
Provides detailed information about estimation results from a haplin object.
}
\usage{
\method{summary}{haplin}(object, reference, ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A \code{haplin} object, i.e. the result of running \code{haplin}.}
  \item{reference}{Same as \code{reference} argument in \code{haplin}. Note that when producing the summary, you can only choose \code{"reciprocal"}, \code{"population"} or \code{"ref.cat"}. You cannot use a numeric value to change the reference category, to do that \code{haplin} must be run over again. (See the \code{reference} argument of \code{haplin}.)}
  \item{...}{Other arguments (ignored).}
}
%\details{}
\references{ Gjessing HK and Lie RT. Case-parent triads: Estimating single- and double-dose effects of fetal and maternal disease gene haplotypes. Annals of Human Genetics (2006) 70, pp. 382-396.\cr\cr
Web Site: \url{http://folk.uib.no/gjessing/genetics/software/haplin/}}
\author{Hakon K. Gjessing\cr
Professor of Biostatistics\cr
Division of Epidemiology\cr
Norwegian Institute of Public Health\cr
\email{hakon.gjessing@uib.no}}
\note{Further information is found on the web page
}
\seealso{\code{\link{haplin}}}

\examples{
\dontrun{

# Produce separate plots for child and mother, dump plots to files:
res <- haplin("data.dat", use.missing = T, maternal = T)
summary(res)

}}




