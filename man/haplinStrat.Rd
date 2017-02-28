\name{haplinStrat}
\alias{haplinStrat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fit \code{haplin} to each subset/stratum of data, determined by the argument \code{strata}}
\description{
 Produces a list, each element of which is an object of class \code{haplin}, which is the result of fitting the log-linear \code{haplin} models to each \code{strata} stratum independently. 
}
\usage{
haplinStrat(filename, data, pedIndex, strata = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{filename}{A character string giving the name and path of the ASCII data file to be read. The file should be in the Haplin data format. }
	\item{data}{An R-object which is the result of using \code{load.gwaa.data} to load data into R. See the web page for a description of how to convert a ped file into a file that can be loaded. The conversion uses \code{prepPed} and \code{convert.snp.ped}.}
	\item{pedIndex}{A file of family indexes constructed by using \code{prepPed} on the original ped file. This file is used by Haplin to extract and store family information.}
	\item{strata}{A single numeric value specifying which data column contains the stratification variable.}
	\item{...}{Remaining arguments to be used by \code{\link{haplin}} in each run.}
}
\details{
	\code{haplinStrat} runs \code{haplin} first on the entire input data file, then on each stratum separately. Strata are defined by the \code{strata} variable, which can be coded as numerical or character. However, one should use only a moderate number of levels/strata, since \code{haplin} will be run independetly on each, and some strata may otherwise have an insufficient amount of data. 
	Running \code{haplinStrat} is thus just a simplification of manually splitting the file into strata and running \code{haplin} on each, with one important difference, however: Since some strata may be small etc., haplin might conceivably choose different haplotypes in different strata, and also choose different reference haplotypes. When first running haplin on the entire file, \code{haplinStrat} saves the selected haplotypes and chosen reference category. Then, in the strata-specific runs haplinStrat forces haplin to choose the same haplotypes/reference category in all runs, so that results from different strata are comparable.
	The main reason for running haplinStrat is to test for gene-environment interactions. This is achieved by running \code{gxe} on the output from \code{haplinStrat}.
	\code{haplinStrat} can also be run from within \code{haplinSlide} by using the \code{strata} argument in \code{haplinSlide}.
}
\value{
	A list of objects of class haplin is returned. The first element contains the result of running \code{haplin} on the entire data file; the remaining elements are the results from each of the strata. The names of the list correspond to the values of the \code{strata} variable.
}
\references{ Gjessing HK and Lie RT. Case-parent triads: Estimating single- and double-dose effects of fetal and maternal disease gene haplotypes. Annals of Human Genetics (2006) 70, pp. 382-396.\cr\cr
Web Site: \url{http://folk.uib.no/gjessing/genetics/software/haplin/}}
\author{Hakon K. Gjessing\cr
Professor of Biostatistics\cr
Division of Epidemiology\cr
Norwegian Institute of Public Health\cr
\email{hakon.gjessing@fhi.no}}

\note{
	Further information is found on the web page.
}
%\seealso{\code{\link{haplin}}, \code{\link{suest}}, \code{\link{summary.haplin}}, \code{\link{plot.haplin}}}
\seealso{\code{\link{haplin}}, \code{\link{summary.haplin}}, \code{\link{plot.haplin}}, \code{\link{haptable}}, \code{\link{toDataFrame}}, \code{\link{haplinSlide}}, \code{\link{gxe}}}
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
