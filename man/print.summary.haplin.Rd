\name{print.summary.haplin}
\alias{print.summary.haplin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Print the summary of a haplin object}
\description{
 Print the result of applying \code{summary} to a haplin object
}
\usage{
\method{print}{summary.haplin}(x, digits, ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A \code{haplin} object, i.e. the result of running \code{haplin}.}
  \item{digits}{The number of digits to be used in the printout. Defaults to 3.}
  \item{...}{Other arguments (ignored).}
}
%\details{}
\references{ Gjessing HK and Lie RT. Case-parent triads: Estimating single- and double-dose effects of fetal and maternal disease gene haplotypes. Annals of Human Genetics (2006) 70, pp. 382-396.\cr\cr
Web Site: \url{http://www.uib.no/smis/gjessing/genetics/software/haplin/}}
\author{Hakon K. Gjessing\cr
Professor of Biostatistics\cr
Division of Epidemiology\cr
Norwegian Institute of Public Health\cr
\email{hakon.gjessing@fhi.no}}
\note{Further information is found on the web page
}
\seealso{\code{\link{haplin}}}

\examples{
\dontrun{

# Standard summary:
res <- haplin("data.dat", use.missing = T, maternal = T)
summary(res)

# Increase number of digits in printout
print(summary(res), digits = 8)

}}





