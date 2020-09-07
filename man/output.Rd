\name{output}
\alias{output}
%\alias{plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Save files with summary, table, and plot from a haplin object.}
\description{
	Create summary tables and figure from a \code{haplin} object. Save results as separate files in a specified directory.
}
\usage{
	output(object, dirname, ask = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A \code{haplin} object, i.e. the result of running \code{haplin}.}
  \item{dirname}{Text string, for instance "c:/work/haplinresults". Name of directory where results should be saved. Default is to save results in the current working directory.}
  \item{ask}{Logical. If TRUE (default), you will be asked before overwriting any files with the same name. If FALSE, \code{output} will overwrite without warning.}
}
\details{ 
	After having run \code{haplin} and saved the result (in the R workspace), the \code{output} function will extract summary results, a summary table, and a plot of the results and save them to the specified directory. The filenames will be haplin_summary.txt, haplin_table.txt and haplin_plot.jpg, respectively. \code{output} simply uses the available functions \code{summary}, \code{haptable}, and \code{plot} to produce the files, but is a quick way of saving all the relevant results.
}
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

# Run haplin and save results in separate files 
# in the c:\work\haplinresults directory:
res <- haplin("data.dat", use.missing = T, maternal = T)
output(res, dirname = "c:/work/haplinresults")
}}




