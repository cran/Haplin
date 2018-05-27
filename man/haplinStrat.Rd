\name{haplinStrat}
\alias{haplinStrat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fit \code{haplin} to each subset/stratum of data, determined by the argument \code{strata}}
\description{
 Produces a list, each element of which is an object of class \code{haplin}, which is the result of fitting the log-linear \code{haplin} models to each \code{strata} stratum independently. 
}
\usage{
haplinStrat( data, strata = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{data}{R-object of class "haplin.ready", which is e.g., output from \code{\link{genDataPreprocess}} or \code{\link{genDataLoad}}, and contains covariate and genetic data.}
	\item{strata}{A single character specifying which data column contains the stratification variable.}
	\item{...}{Remaining arguments to be used by \code{\link{haplin}} in each run.}
}
\details{
	\code{haplinStrat} runs \code{haplin} first on the entire input data file, then on each stratum separately. Strata are defined by the \code{strata} variable, which can be coded as numerical or character. However, one should use only a moderate number of levels/strata, since \code{haplin} will be run independetly on each, and some strata may otherwise have an insufficient amount of data. 
	Running \code{haplinStrat} is thus just a simplification of manually splitting the file into strata and running \code{haplin} on each; the end result would be the same.
	The main reason for running haplinStrat is to test for gene-environment interactions. This is achieved by running \code{postTest} on the result from \code{haplinStrat}.
	\code{haplinStrat} can also be run from within \code{haplinSlide} by using the \code{strata} argument in \code{haplinSlide}.
}
\value{
	A list of objects of class haplin is returned. The first element contains the result of running \code{haplin} on the entire data file; the remaining elements are the results from each of the strata. The names of the list correspond to the values of the \code{strata} variable.
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
\seealso{\code{\link{haplin}}, \code{\link{summary.haplin}}, \code{\link{plot.haplin}}, \code{\link{haptable}}, \code{\link{toDataFrame}}, \code{\link{haplinSlide}}, \code{\link{gxe}}}
\examples{

# setting up the directory with exemplary data
dir.in <- system.file( "extdata", package = "Haplin" )
file.in <- paste0( dir.in, "/data.dat" )

# reading data in
data.in <- genDataRead( file.in, file.out = "poo_exmpl_data_read", format = "haplin",
  n.vars = 1, allele.sep = " ", col.sep = " ", overwrite = TRUE )
# preprocessing the data
data.preproc <- genDataPreprocess( data.in, design = "triad",
  file.out = "poo_exmpl_data_preproc", overwrite = TRUE )

# running haplinStrat, checking for gene-environment interactions
res.GxE <- haplinStrat( data.preproc, markers = 2, strata = 1, poo = FALSE,
  response = "mult", reference = 2, use.missing = TRUE )
res.GxE

# running haplinStrat, checking for POO-environment interactions
res.POOxE <- haplinStrat( data.preproc, markers = c(1,2,3), strata = 1, poo = TRUE,
  response = "mult", reference = "ref.cat", use.missing = TRUE )
res.POOxE

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
postTest(result)
%# Compute an overall p-value for the scan, corrected for multiple testing
%# and dependencies between windows:
%suest(result.1)

}

}
