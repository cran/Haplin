\name{haplin}
\alias{haplin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fitting log-linear models to case-parent triad and/or case-control data}
\description{
 Produces an object of class \code{haplin}, which is the result of fitting the log-linear models to the data 
}
\usage{
haplin(filename, 
markers = "ALL", n.vars = 0, sep = " ", allele.sep = ";", 
na.strings = "NA", design = "triad", use.missing = FALSE, 
xchrom = FALSE, maternal = FALSE, test.maternal = FALSE, 
scoretest = "no", ccvar = NULL, covar = NULL, sex = NULL, 
reference = "reciprocal", response = "free", 
threshold = 0.01, max.haplos = NULL, haplo.file = NULL, 
resampling = FALSE, max.EM.iter = 50, data.out = "no", 
verbose = TRUE, printout = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{}{Of the following arguments, only \code{filename} is required. Use of the remaining arguments will depend on the type of analysis.}
	\item{filename}{A character string giving the name and path of the ASCII data file to be read. }
	\item{markers}{ Default is "ALL", which means HAPLIN uses all available markers in the data set in the analysis. For the current version of HAPLIN the number of markers used at a single run should probably not exceed 4 or 5 due to the computational burden. The markers argument can be used to select appropriate markers from the file without creating a new file for the selected markers. For instance, if markers is set to c(2,4), HAPLIN will only use the second and fourth markers supplied in the data set. When running HAPLIN, it may be a good idea to start exploring a few markers at a time, using this argument.  }
	\item{n.vars}{Numeric. The number of variables (columns) in the data file before (to the left) of the genetic data.}
	\item{sep}{The character separator used in the data file to separate between "columns", where each column contains the two alleles of a single individual at a single marker. }
	\item{allele.sep}{ The character separator used in the data file to separate the two alleles for a single individual in a single marker. The recommended (default) separator is ";", but for SNPs an empty "" is also common.}
	\item{na.strings}{ The character string indicating missing data in the data file. Default is to use "NA" in place of, for instance, C;T for a SNP that hasn't been typed in that individual.}
	\item{design}{The value "triad" is used for the standard case triad design, without indepdendent controls. The value "cc.triad" means a combination of case triads and control triads. This requires the argument \code{ccvar} to point to the data column containing the case-control variable. The value "cc" means a simple case-control design, where the parents have not been genotyped (there are no data columns for parental genes).}
	\item{use.missing}{ A logical value used to determine whether triads with missing data should be included in the analysis. When set to TRUE, Haplin uses the EM algorithm to obtain risk estimates, also taking into account triads with missing data. The standard errors and p-values are adjusted to correct for this. The default, however, is FALSE. When FALSE, all triads having  any sort of missing data are excluded before the analysis is run. Note that Haplin only looks at markers actually used in the analysis,  so that if the markers  argument (see below) is used to select a collection of markers for analysis, Haplin only excludes triads with missing data on the included markers. }
	\item{xchrom}{Logical, defaults to "FALSE". If set to "TRUE", haplin assumes the markers are on the x-chromosome. This option should be combined with specifying the \code{sex} argument, and setting (for the time being) \code{response = "mult", reference = "ref.cat", maternal = F}.}
	\item{maternal}{ If TRUE, maternal effects are estimated as well as the standard fetal effects. }
	\item{test.maternal}{Not yet implemented.}
	\item{scoretest}{Special interest only. If "no", no score test is computed. If "yes", an overall score p-value is included in the output, and the individual score values are returned in the haplin object. If "only", haplin is only run under the null hypothesis, and a simple score object is returned instead of the full haplin object. Useful if only score testing is needed.}
	\item{ccvar}{Numeric. Should give the column number for the column containing the case-control indicator in the data file. Needed for the "cc" and "cc.triad" designs. The column should contain two numeric values, of which the largest one is always used to denote cases.}
	\item{covar}{Not yet implemented.}
	\item{sex}{A numeric value specifying which of the data columns that contains the sex variable. The variable should be coded 1 for males and 2 for females. To be used with \code{xchrom = TRUE}.}
	\item{reference}{ Decides how HAPLIN chooses its reference category for the effect estimates. Default value is "reciprocal". With the reciprocal reference the effect of a single or double dose of each haplotype is measured relative to the remaining haplotypes. This means that a new reference category is used for each single haplotype. Other possible values are "population" (which is similar to reciprocal, but where the reference category is always the total population), and "ref.cat", where a single haplotype is used as reference for all the rest. For ref.cat, the default is to choose the most frequent haplotype as the reference haplotype. The reference haplotype can be set explicitly by giving a numeric value for the reference argument. Note that the numeric value refers to the haplotype's position among the haplotypes selected for analysis by HAPLIN. This means that one should run HAPLIN once first to see what haplotypes are used before giving a numeric value to reference. }
	\item{response}{The default value "free" means that both single- and double dose effects are estimated. Choosing "mult" instead specifies a multiplicative dose-response model. }
	\item{threshold}{ Sets the (approximate) lower limit for the haplotype frequencies of those haplotypes that should be retained in the analysis. Hapotypes that are less frequent are removed, and information about this is given in the output. }
	\item{max.haplos}{Not yet implemented.}
	\item{haplo.file}{Not yet implemented.}
	\item{resampling}{ Default is FALSE. When FALSE, the individual haplotypes reconstructed by the EM algorithm as assumed known when computing CIs and p-values. If set to "jackknife" a jackknife-based resampling procedure is used when computing confidence intervals and p-values for effect estimates. This takes more time, but corrects the CIs and p-values for the uncertainty contained in unphased data. Note: in this version of Haplin, the resampling is no longer needed since the confidence intervals and p-values are already corrected in the standard computation. }
	\item{max.EM.iter}{ The maximum number of iterations used by the EM algorithm. This value can be increased if necessary, which sometimes is the case with e.g. case-control data which a substantial amount of missing. However, for triad data with little missing information there is usually no need for many iterations. }
	\item{data.out}{Character. Accepts values "no", "prelim", "null" or "full", with "no" as default. For values other than default, \code{haplin} returns the data file prepared for analysis rather than the usual \code{haplin} estimation results. The data file contains the haplotypes identified for each triad, and a vector of weights giving the probability distribution of different haplotype configurations within a triad. The probabilities are computed from preliminary haplotype frequency estimates, from the null model or from the full likelihood model. The "prelim" option will be much faster but somewhat less precise than the likelihood models.}
	\item{verbose}{ Default is T (=TRUE). During the EM algorithm, HAPLIN prints the estimated parameters and deviance for each step. To avoid the output, set this argument to F (=FALSE). }
	\item{printout}{Logical. If TRUE (default), haplin prints a full summary of the results after finishing the estimation. If FALSE, no such printout is given, but the \code{summary} function can later be applied to a saved result to get the same summary.}
}
\details{
  The output can be examined by print, summary and plot.
}
\value{
  An object of class \code{haplin} is returned.
 (The only exception is when \code{data.out} is set different from "no", where \code{haplin} will produce a data file with haplotypes identified.)
}
\references{ Gjessing HK and Lie RT. Case-parent triads: Estimating single- and double-dose effects of fetal and maternal disease gene haplotypes. Annals of Human Genetics (2006) 70, pp. 382-396.\cr\cr
Web Site: \url{http://www.uib.no/smis/gjessing/genetics/software/haplin/}}
\author{Hakon K. Gjessing\cr
Professor of Biostatistics\cr
Division of Epidemiology\cr
Norwegian Institute of Public Health\cr
\email{hakon.gjessing@fhi.no}}

\section{Warning}{Typically, some of the included haplotypes will be relatively rare, such as a frequency of 1\% - 5\%. For those haplotypes there may be too little data to estimate the double doses properly, so the estimates may be unreliable. This is seen from the extremely wide confidence intervals. The rare double dose estimates should be disregarded, but the remaining single and double dose estimates are valid. To avoid the problem one can also reduce the model to a purely multiplicative model by setting \code{response = "mult"} combined with \code{reference = "ref.cat"}.}

\note{Further information is found on the web page.
}
\seealso{\code{\link{summary.haplin}}, \code{\link{plot.haplin}}, \code{\link{pedToHaplin}}, \code{\link{haplinSlide}}}
\examples{

\dontrun{

# Standard run:
haplin("data.dat")

# Specify path, estimate maternal effects:
haplin("C:/work/data.dat", maternal = T)

# Specify path, use haplotype no. 2 as reference:
haplin("C:/work/data.dat", reference = 2)

# Remove more haplotypes from estimation by increasing the threshold 
# to 5\%:
haplin("C:/work/data.dat", threshold = 0.05)

# Estimate maternal effects, using the most frequent haplotype as reference. 
# Use all data, including triads with missing data. Select 
# markers 3, 4 and 8 from the supplied data.
haplin("C:/work/data.dat", use.missing = T, maternal = T, 
reference = "ref.cat", markers = c(3,4,8))
# Note: in this version of Haplin, the jackknife is 
# no longer necessary since the standard errors are already corrected.

# Some examples showing how to save the Haplin result and later 
# recall plot and summary results:

# Same analysis as above, saving the result in the object "result.1":
result.1 <- haplin("C:/work/data.dat", use.missing = T, maternal = T, 
reference = "ref.cat", markers = c(3,4,8))

# Replot the saved result (fetal effects):
plot(result.1)

# Replot the saved result (maternal effects):
plot(result.1, plot.maternal = T)

# Print a very short summary of saved result:
result.1

# A full summary of saved result, with confidence intervals and 
# p-values (the same as haplin prints when running):
summary(result.1)

# Some examples when the data file contains two covariates, 
# the second is the case-control variable:

# The following standard triad run is INCORRECT since it disregards 
# case status:
haplin("data.dat", use.missing = T, n.vars = 2, design = "triad")

# Combined run on "hybrid" design, correctly using both case-parent 
# triads and control-parent triads:
haplin("data.dat", use.missing = T, n.vars = 2, ccvar = 2, 
design = "cc.triad")

# If parent columns are not in the file, a plain case-control 
# run can be used:
haplin("data.dat", use.missing = T, n.vars = 2, ccvar = 2, 
design = "cc", response = "mult", reference = "ref.cat")

# An example of how to produce a data file with all possible haplotypes
# identified for each triad, together with their probaility weights:
result.data <- haplin("C:/work/data.dat", use.missing = T, 
markers = c(3,4,8), data.out = "prelim")
# result.data will then contain the data file, with a vector of 
# probabilities (freq) computed from the preliminary haplotype
# frequencies.



}

}
