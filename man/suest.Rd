\name{suest}
\alias{suest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute a joint p-value for a list of haplin fits (usually from a sliding window approach), correcting for multiple testing.}
\description{
 The first argument to \code{suest} should be a list of haplin estimation results (from the same data file), usually the output from \code{\link{haplinSlide}}. \code{suest} produces as a result a joint overall p-value based on aggregating the individual p-values and then correcting for multiple testing. The correction is achieved by using the principle of "seemingly unrelated" estimation, taking into account the correlation between the individual estimation results. 
}
\usage{
suest(reslist)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{reslist}{A list whose elements are different haplin runs on the same data file, typically the output of \code{\link{haplinSlide}}.}
}
\details{
	\code{\link{haplinSlide}} runs \code{haplin} on a series of overlapping
    windows of markers from the same data file, typically within the same gene.
    Since each run produces a separate overall p-value, \code{suest} computes a
    joint overall p-value for the gene (or region) that has been scanned. It
    corrects the overall p-value for multiple testing, also taking into account
    the fact that the sequence of estimates produced by \code{haplinSlide} will
    be dependent, both because they are computed on the same data set and also
    since the windows are overlapping (if the window length is larger than 1).
    If the \code{suest} estimation fails (which doesn't happen very often), a
    standard Bonferroni correction is used instead. Important:
    \code{haplinSlide} must be run with the option \code{table.output = FALSE}
    to provide \code{suest} with enough information.
}
\value{
	A list is returned, the most important elements of which are:
	\item{pval.obs}{The overall score p-values from each haplin run}
	\item{pval.obs.corr}{The joint p-value, corrected for multiple testing}
	\item{bonferroni}{A logical, usually FALSE, which means the \code{suest} estimation went well. If TRUE, it means that the \code{suest} estimation failed for some reason, and a standard Bonferroni correction was used instead.}
}
\references{ Gjessing HK and Lie RT. Case-parent triads: Estimating single- and double-dose effects of fetal and maternal disease gene haplotypes. Annals of Human Genetics (2006) 70, pp. 382-396.\cr\cr
Web Site: \url{http://folk.uib.no/gjessing/genetics/software/haplin/}}
\author{Hakon K. Gjessing\cr
Professor of Biostatistics\cr
Division of Epidemiology\cr
Norwegian Institute of Public Health\cr
\email{hakon.gjessing@uib.no}}

\note{
	Further information is found on the web page.
}
\seealso{\code{\link{haplin}}, \code{\link{haplinSlide}}}
\examples{

\dontrun{
# (Almost) all standard haplin runs can be done with haplinSlide. 
# Below is an illustration. See the haplin help page for more 
# examples.
# 
# Analyzing the effect of fetal genes, including triads with missing data,
# using a multiplicative response model. When winlength = 1, separate
# markers are used. To make longer windows, winlength can be increased
# correspondingly:
result.1 <- haplinSlide("C:/work/data.dat", use.missing = T, response = "mult",
reference = "ref.cat", winlength = 1, table.output = F)
# Provide summary of separate results:
lapply(result.1, summary)
# Plot results:
par(ask = T)
lapply(result.1, plot)
# Compute an overall p-value for the scan, corrected for multiple testing
# and dependencies between windows:
suest(result.1)

}

}
