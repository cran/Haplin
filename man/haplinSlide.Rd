\name{haplinSlide}
\alias{haplinSlide}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Run haplin analysis in a series of sliding windows over a sequence of markers/SNPs}
\description{
 Produces a list, each element of which is an object of class \code{haplin}, which is the result of fitting the log-linear \code{haplin} models to the data one "window" at a time. 
}
\usage{
haplinSlide( data, markers = "ALL", winlength = 1, 
strata = NULL, table.output = TRUE, cpus = 1, para.env = NULL, slaveOutfile = "", 
printout = FALSE, verbose = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{data}{R-object of class "haplin.ready", which is e.g., output from \code{\link{genDataPreprocess}} or \code{\link{genDataLoad}}, and contains covariate and genetic data.}
	\item{markers}{Default is "ALL", which means haplinSlide uses all available markers in the data set in the analysis. Alternatively, the relevant markers can be specified by, for instance, \code{markers = c(1, 3:10)}, which would use the 10 first markers except marker 2. \code{haplinSlide} will then run haplin on a series of windows selected from the supplied \code{markers}. The \code{winlength} argument decides the length of the windows. See details.}
	\item{winlength}{Length of the sliding, overlapping windows to be run along the markers. See details.}
	\item{strata}{A single numeric value specifying which data column contains the stratification variable.}
	\item{table.output}{If \code{TRUE}, the \code{haptable} function will be applied to each result after estimation, greatly reducing the size of the output. If \code{FALSE}, each element of the output list is a standard \code{haplin} object. To preserve memory, default is set to TRUE.}
	\item{cpus}{\code{haplinSlide} allows parallel processing of its analyses. The \code{cpus} argument should preferably be set to the number of available cpu's. If set lower, it will save some capacity for other processes to run. Setting it too high should not cause any serious problems.}
	\item{para.env}{The user can choose parallel environment to use --- "parallel" (default) or "Rmpi" (for use on clusters); this option is used only when \code{cpus} argument is larger than 1.}
	\item{slaveOutfile}{Character. To be used when \code{cpus > 1}. If \code{slaveOutfile = ""} (default), output from all running cores will be printed in the standard R session window. Alternatively, the output can be saved to a file by specifying the file path and name.}
	\item{printout}{Default is FALSE. If TRUE, provides a full summary of each \code{haplin} result during the run of \code{haplinSlide}.}
	\item{verbose}{Same as for \code{haplin}, but defaults to FALSE to reduce output size.}
	\item{...}{Remaining arguments to be used by \code{\link{haplin}} in each run.}
}
\details{
	\code{haplinSlide} runs \code{haplin} on a series of overlapping windows of the chosen markers. Except for the \code{markers} and \code{winlength} arguments, all arguments are used exactly as in \code{haplin} itself. For instance, if \code{markers = c(1, 3, 4, 5, 7, 8)} and \code{winlength = 4}, \code{haplinSlide} will run \code{haplin} on first the markers \code{c(1, 3, 4, 5)}, then on \code{c(3, 4, 5, 7)}, and finally on \code{c(4, 5, 7, 8)}. The results are returned in a list. The elements are named "1-3-4-5" etc., and can be extracted with, say, \code{summary(res[["1-3-4-5"]])} etc., where \code{res} is the saved result. Or the output can be examined by, for instance, using \code{lapply(res, summary)} and \code{lapply(res, plot)}.\cr
	When running \code{haplinSlide} on a large number of markers, the output can become prohibitively large. In that case \code{table.output} should be set to \code{TRUE}, and \code{haplinSlide} will return a list of summary "haptables". This list can then be stacked into a single dataframe using \code{toDataFrame}. To avoid exessive memory use, the default is \code{table.output = TRUE}.\cr
	When multiple cores are available, set the \code{cpus} to the number of cores that should be used. This will run \code{haplinSlide} in parallel on the chosen number of cores. Note that feedback is provided by each of the cores separately, and some cores may start working on markers far out in the sequence.
}
\value{
	A list of objects of class haplin is returned.
}
\references{ Gjessing HK and Lie RT. Case-parent triads: Estimating single- and double-dose effects of fetal and maternal disease gene haplotypes. Annals of Human Genetics (2006) 70, pp. 382-396.\cr\cr
Web Site: \url{http://www.uib.no/smis/gjessing/genetics/software/haplin/}}
\author{Hakon K. Gjessing\cr
Professor of Biostatistics\cr
Division of Epidemiology\cr
Norwegian Institute of Public Health\cr
\email{hakon.gjessing@uib.no}}

\note{
	Further information is found on the web page.
}
%\seealso{\code{\link{haplin}}, \code{\link{suest}}, \code{\link{summary.haplin}}, \code{\link{plot.haplin}}}
\seealso{\code{\link{haplin}}, \code{\link{summary.haplin}}, \code{\link{plot.haplin}}, \code{\link{haptable}}, \code{\link{toDataFrame}}}
\examples{

\dontrun{
# (Almost) all standard haplin runs can be done with haplinSlide. 
# Below is an illustration. See the haplin help page for more 
# examples.
# 

# 1. Read the data:
my.haplin.data <- genDataRead( file.in = "HAPLIN.trialdata.txt", file.out =
  "trial_data1", dir.out = ".", format = "haplin", n.vars = 0 )

# 2. Run pre-processing:
haplin.data.prep <- genDataPreprocess( data.in = my.haplin.data,
  format = "haplin", design = "triad", file.out = "trial_data1_prep",
  dir.out = "." )

# 3. Analyze:
# Analyzing the effect of fetal genes, including triads with missing data,
# using a multiplicative response model. When winlength = 1, separate
# markers are used. To make longer windows, winlength can be increased
# correspondingly:
result.1 <- haplinSlide( haplin.data.prep, use.missing = T, response = "mult",
reference = "ref.cat", winlength = 1, table.output = F)
# Provide summary of separate results:
lapply(result.1, summary)
# Plot results:
par(ask = T)
lapply(result.1, plot)
%# Compute an overall p-value for the scan, corrected for multiple testing
%# and dependencies between windows:
%suest(result.1)

}

}
