\name{gxe}
\alias{gxe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Test for gene-environment interaction}
\description{
 Performes a gene-environment test to check if \code{haplin} estimates of relative risks change over strata of environment. It is typically applied to the output from \code{haplinStrat}
}
\usage{
gxe(object.list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{object.list}{A list of \code{haplin} results, almost always the output from \code{haplinStrat}. The first element is the result of running \code{haplin} on all data; the remaining elements are the results for each stratum separately.}
}
\details{
	\code{haplinStrat} runs \code{haplin} first on the entire input data file, then on each stratum separately. The results from \code{haplinStrat} are similar to just manually splitting the file into strata and running \code{haplin} on each, with one important difference, however: Since some strata may be small etc., haplin might conceivably choose different haplotypes in different strata, and also choose different reference haplotypes. When first running haplin on the entire file, \code{haplinStrat} saves the selected haplotypes and chosen reference category. Then, in the strata-specific runs haplinStrat forces haplin to choose the same haplotypes/reference category in all runs, so that results from different strata are comparable.
	When applying \code{gxe} to the output from haplinStrat, it will test whether there is a statistically significant change in parameter estimates from stratum to stratum, i.e. a gene-environment interaction since strata usually are defined by an environmental exposure.
	\code{gxe} uses Wald tests to test for interactions. It always tests whether there is change in haplotype frequencies from stratum to stratum. More importantly, it separately tests whether any genetic effects, such as fetal genetic effects, maternal effects, or parent-of-origin effects, change significantly over strata.
	\code{gxe} can also be run from within \code{haplinSlide} by using the \code{strata} argument in \code{haplinSlide}.
}
\value{
	A dataframe with one row for each test that is performed (haplo.freq is the first, the remaining depend on the model that has been estimated). The Wald chi-squared test value, degrees-of-freedom, and resulting p-value are reported.
}
\references{ Gjessing HK and Lie RT. Case-parent triads: Estimating single- and double-dose effects of fetal and maternal disease gene haplotypes. Annals of Human Genetics (2006) 70, pp. 382-396.\cr\cr
Web Site: \url{https://people.uib.no/gjessing/genetics/software/haplin/}}
\author{Hakon K. Gjessing\cr
Professor of Biostatistics\cr
Division of Epidemiology\cr
Norwegian Institute of Public Health\cr
\email{hakon.gjessing@uib.no}}
\section{NOTE:}{In the future, the structure of the output from gxe will change. In particular, measures of ratios of relative risks will be reported in addition to the p-values}
\note{
	Further information is found on the web page.
}
%\seealso{\code{\link{haplin}}, \code{\link{suest}}, \code{\link{summary.haplin}}, \code{\link{plot.haplin}}}
\seealso{\code{\link{haplin}}, \code{\link{haplinStrat}}, \code{\link{haplinSlide}}}
\examples{

\dontrun{
# All standard haplin runs can be done with haplinStrat. 
# Below is an illustration. See the haplin help page for more 
# examples.
# 
# Analyzing the effect of fetal genes, including triads with missing data,
# using a multiplicative response model. The first column of the data file
# in this example contains the stratification variable.
result <- haplinStrat("C:/work/data.dat", strata = 1, use.missing = T, response = "mult",
reference = "ref.cat", winlength = 1)
# Provide summary of separate results:
lapply(result, summary)
# Plot results separately:
par(ask = T)
lapply(result, plot)
#
# Convert results to table format and stack them over strata:
haptable(result)
# Test for interaction between haplotype risk estimates and the strata variable:
gxe(result)
%# Compute an overall p-value for the scan, corrected for multiple testing
%# and dependencies between windows:
%suest(result.1)

}

}
