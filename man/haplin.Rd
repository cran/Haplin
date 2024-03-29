\name{haplin}
\alias{haplin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fitting log-linear models to case-parent triad and/or case-control data}
\description{
\code{haplin} fits a log-linear model to case-parent triads, case-control data, or combined (hybrid) case-parent control-parent triads or dyads. It estimates marker or haplotype frequencies, and uses the EM algorithm to reconstruct haplotypes and, if requested, impute missing genotypes. \code{haplin} prints and plots estimates of relative risks associated with fetal and maternal haploypes, and in addition allows splitting fetal haplotype effects into maternally and paternally inherited effects. It allows special models, like x-inactivation, to be fitted on the X-chromosome. The result is an object of class \code{haplin}, which can be explored with \code{summary}, \code{plot}, and \code{haptable}.
}
\usage{
haplin( data, markers = "ALL", 
design = "triad", use.missing = FALSE, 
xchrom = FALSE, maternal = FALSE, test.maternal = FALSE, 
poo = FALSE, scoretest = "no", ccvar = NULL, strata = NULL, 
sex = NULL, comb.sex = "double",
reference = "reciprocal", response = "free", 
threshold = 0.01, max.haplos = NULL, haplo.file = NULL, 
resampling = "no", max.EM.iter = 50, data.out = "no", 
verbose = TRUE, printout = TRUE )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{data}{An R-object which is the result of using \link{genDataPreprocess}. See the web page for a detailed description of how to use this function.}
	\item{markers}{Default is "ALL", which means \code{haplin} uses all available markers in the data set in the analysis. For the current version of \code{haplin} the number of markers used at a single run should probably not exceed 4 or 5 due to the computational burden. The markers argument can be used to select appropriate markers from the file without creating a new file for the selected markers. The relevant markers can be specified by giving a vector or numbers (e.g., \code{markers = c(1, 3:10)} will use the 10 first markers except marker 2) or characters (e.g., \code{markers = c("m1", "m3", "rs35971")}). When running \code{haplin}, it may be a good idea to start exploring a few markers at a time, using this argument.  }
	\item{design}{The value "triad" is used for the standard case triad design, without independent controls. The value "cc.triad" means a combination of case triads and control triads. This requires the argument \code{ccvar} to point to the data column containing the case-control variable. The value "cc" means a simple case-control design, where the parents have not been genotyped (there are no data columns for parental genes). NOTE: \code{design} is also set in \code{genDataPreprocess}. Almost always, the two arguments should be equal. Occasionally, however, the user might want to override the original argument by switching from 'cc.triad' to 'triad' or vice versa.	}
	\item{use.missing}{ A logical value used to determine whether triads with missing data should be included in the analysis. When set to TRUE, \code{haplin} uses the EM algorithm to obtain risk estimates, also taking into account triads with missing data. The standard errors and p-values are adjusted to correct for this. The default, however, is FALSE. When FALSE, all triads having  any sort of missing data are excluded before the analysis is run. Note that \code{haplin} only looks at markers actually used in the analysis,  so that if the markers  argument (see below) is used to select a collection of markers for analysis, \code{haplin} only excludes triads with missing data on the included markers. }
	\item{xchrom}{Logical, defaults to "FALSE". If set to "TRUE", \code{haplin} assumes the markers are on the x-chromosome. This option should be combined with specifying the \code{sex} argument. In addition, \code{comb.sex} can be useful. \code{xchrom = T} can be combined with \code{poo = T} and/or \code{maternal = T}.}
	\item{maternal}{ If TRUE, maternal effects are estimated as well as the standard fetal effects. }
	\item{test.maternal}{Not yet implemented.}
	\item{poo}{Parent-Of-Origin effects. If TRUE, \code{haplin} will split single-dose effects into two separate effect estimates, one for the maternally inherited haplotype, and one for the paternally inherited haplotype. Double dose will be estimated as before.}
	\item{scoretest}{Special interest only. If "no", no score test is computed. If "yes", an overall score p-value is included in the output, and the individual score values are returned in the \code{haplin} object. If "only", \code{haplin} is only run under the null hypothesis, and a simple score object is returned instead of the full \code{haplin} object. Useful if only score testing is needed.}
	\item{ccvar}{Numeric. Should give the column number for the column containing the case-control indicator in the data file. Needed for the "cc" and "cc.triad" designs. The column should contain two numeric values, of which the largest one is always used to denote cases.}
	\item{strata}{Not yet implemented.}
	\item{sex}{To be used with \code{xchrom = TRUE}. A numeric value specifying which of the data columns that contains the sex variable. The variable should be coded 1 for males and 2 for females.}
	\item{comb.sex}{To be used with \code{xchrom = TRUE}. A character value that specifies how to handle gender differences on the X-chromosome. If set to "males" or "females", analyses are done either for just males or just females, respectively. If set to "single" or "double", males and females are used in a combined analysis. Specifically, when "single", the effect of a (single) allele in males is assumed to equal the effect of a single allele dose in females, and similarly, when "double", a single allele in males is assumed to have the same effect as a double allele dose in females. Default is "double", which corresponds to X-inactivation. See separate description for more details.}
	\item{reference}{Decides how \code{haplin} chooses its reference category for the effect estimates. Default value is "reciprocal". With the reciprocal reference the effect of a single or double dose of each haplotype is measured relative to the remaining haplotypes. This means that a new reference category is used for each single haplotype. Other possible values are "population" (which is similar to reciprocal, but where the reference category is always the total population), and "ref.cat", where a single haplotype is used as reference for all the rest. For ref.cat, the default is to choose the most frequent haplotype as the reference haplotype. The reference haplotype can be set explicitly by giving a numeric value for the reference argument. Note that the numeric value refers to the haplotype's position among the haplotypes selected for analysis by \code{haplin}. This means that one should run \code{haplin} once first to see what haplotypes are used before giving a numeric value to reference. }
	\item{response}{The default value "free" means that both single- and double dose effects are estimated. Choosing "mult" instead specifies a multiplicative dose-response model. }
	\item{threshold}{Sets the (approximate) lower limit for the haplotype frequencies of those haplotypes that should be retained in the analysis. Hapotypes that are less frequent are removed, and information about this is given in the output. Default is 0.01.}
	\item{max.haplos}{Not yet implemented.}
	\item{haplo.file}{Not yet implemented.}
	\item{resampling}{Mostly for testing. Default is "no". When "no", the individual haplotypes reconstructed by the EM algorithm as assumed known when computing CIs and p-values. If set to "jackknife" a jackknife-based resampling procedure is used when computing confidence intervals and p-values for effect estimates. This takes more time, but corrects the CIs and p-values for the uncertainty contained in unphased data. Note: in all recent versions of \code{haplin}, the resampling is no longer needed since the confidence intervals and p-values are already corrected in the standard computation. }
	\item{max.EM.iter}{The maximum number of iterations used by the EM algorithm. This value can be increased if necessary, which sometimes is the case with e.g. case-control data which a substantial amount of missing. However, for triad data with little missing information there is usually no need for many iterations. }
	\item{data.out}{Character. Accepts values "no", "prelim", "null" or "full", with "no" as default. For values other than default, \code{haplin} returns the data file prepared for analysis rather than the usual \code{haplin} estimation results. The data file contains the haplotypes identified for each triad, and a vector of weights giving the probability distribution of different haplotype configurations within a triad. The probabilities are computed from preliminary haplotype frequency estimates, from the null model or from the full likelihood model. The "prelim" option will be much faster but somewhat less precise than the likelihood models.}
	\item{verbose}{ Default is T (=TRUE). During the EM algorithm, \code{haplin} prints the estimated parameters and deviance for each step. To avoid the output, set this argument to F (=FALSE). }
	\item{printout}{Logical. If TRUE (default), \code{haplin} prints a full summary of the results after finishing the estimation. If FALSE, no such printout is given, but the \code{summary} function can later be applied to a saved result to get the same summary.}
}
\details{
Input data can be either a \code{haplin} format data file, or a PED data. These have to be loaded into R first, using \link{genDataRead} or \link{genDataLoad} functions, and then pre-processed with the \link{genDataPreprocess} function. If the PED data file is used, the arguments \code{filename}, \code{n.vars}, \code{sep}, \code{allele.sep}, \code{na.strings}, \code{ccvar}, and \code{sex} need not be specified.\cr
The output can be examined by \code{print}, \code{summary}, \code{plot} and \code{haptable}.
}
\value{
  An object of class \code{haplin} is returned.
 (The only exception is when \code{data.out} is set different from "no", where \code{haplin} will produce a data file with haplotypes identified.)
}
\references{ Gjessing HK and Lie RT. Case-parent triads: Estimating single- and double-dose effects of fetal and maternal disease gene haplotypes. Annals of Human Genetics (2006) 70, pp. 382-396.\cr\cr
Web Site: \url{https://haplin.bitbucket.io}}
\author{Hakon K. Gjessing\cr
Professor of Biostatistics\cr
Division of Epidemiology\cr
Norwegian Institute of Public Health\cr
\email{hakon.gjessing@uib.no}}

\section{Warning}{Typically, some of the included haplotypes will be relatively rare, such as a frequency of 1\% - 5\%. For those haplotypes there may be too little data to estimate the double doses properly, so the estimates may be unreliable. This is seen from the extremely wide confidence intervals. The rare double dose estimates should be disregarded, but the remaining single and double dose estimates are valid. To avoid the problem one can also reduce the model to a purely multiplicative model by setting \code{response = "mult"} combined with \code{reference = "ref.cat"}.}

\note{Further information is found on the web page.
}
\seealso{\code{\link{summary.haplin}}, \code{\link{plot.haplin}}, \code{\link{pedToHaplin}}, \code{\link{haptable}}, \code{\link{haplinSlide}}, \code{\link{genDataLoad}}, \code{\link{genDataRead}}, \code{\link{genDataPreprocess}} }
\examples{
# setting up the directory with exemplary data
dir.in <- system.file( "extdata", package = "Haplin" )
file.in <- file.path( dir.in, "data.dat" )

# reading data in
data.in <- genDataRead( file.in, file.out = "poo_exmpl_data_read", format = "haplin",
  dir.out = tempdir( check = TRUE ), n.vars = 1, allele.sep = " ", col.sep = " ",
  overwrite = TRUE )
# preprocessing the data
data.preproc <- genDataPreprocess( data.in, design = "triad",
  file.out = "poo_exmpl_data_preproc", dir.out = tempdir( check = TRUE ), overwrite = TRUE )

# running haplin, calculating POO
res.POO <- haplin( data.preproc, markers = 2, poo = TRUE, response = "mult",
  reference = 2, use.missing = TRUE )
res.POO

\dontrun{
# 1. Read the data:
my.haplin.data <- genDataRead( file.in = "HAPLIN.trialdata.txt", file.out =
  "trial_data1", dir.out = ".", format = "haplin", n.vars = 0 )

# 2. Run pre-processing:
haplin.data.prep <- genDataPreprocess( data.in = my.haplin.data, format =
  "haplin", design = "triad", file.out = "trial_data1_prep", dir.out = "." )

# 3. Analyze:
# Standard run:
haplin( haplin.data.prep )

# Specify path, estimate maternal effects:
haplin( haplin.data.prep, maternal = T )

# Specify path, use haplotype no. 2 as reference:
haplin( haplin.data.prep, reference = 2 )

# Remove more haplotypes from estimation by increasing the threshold 
# to 5\%:
haplin( haplin.data.prep, threshold = 0.05 )

# Estimate maternal effects, using the most frequent haplotype as reference. 
# Use all data, including triads with missing data. Select 
# markers 3, 4 and 8 from the supplied data.
haplin( haplin.data.prep, use.missing = T, maternal = T, 
reference = "ref.cat", markers = c(3,4,8) )
# Note: in this version of haplin, the jackknife is 
# no longer necessary since the standard errors are already corrected.

# Some examples showing how to save the haplin result and later 
# recall plot and summary results:

# Same analysis as above, saving the result in the object "result.1":
result.1 <- haplin( haplin.data.prep, use.missing = T, maternal = T, 
reference = "ref.cat", markers = c(3,4,8) )

# Replot the saved result (fetal effects):
plot( result.1 )

# Replot the saved result (maternal effects):
plot( result.1, plot.maternal = T )

# Print a very short summary of saved result:
result.1

# A full summary of saved result, with confidence intervals and 
# p-values (the same as haplin prints when running):
summary( result.1 )

# Some examples when the data file contains two covariates, 
# the second is the case-control variable:

# The following standard triad run is INCORRECT since it disregards 
# case status:
haplin("data.dat", use.missing = T, n.vars = 2, design = "triad")

# Combined run on "hybrid" design, correctly using both case-parent 
# triads and control-parent triads:
haplin( my.haplin.data, use.missing = T, n.vars = 2, ccvar = 2, 
design = "cc.triad" )

# If parent columns are not in the file, a plain case-control 
# run can be used:
haplin( my.haplin.data, use.missing = T, n.vars = 2, ccvar = 2, 
design = "cc", response = "mult", reference = "ref.cat" )

# An example of how to produce a data file with all possible haplotypes
# identified for each triad, together with their probaility weights:
result.data <- haplin( my.haplin.data, use.missing = T, 
markers = c(3,4,8), data.out = "prelim" )
# result.data will then contain the data file, with a vector of 
# probabilities (freq) computed from the preliminary haplotype
# frequencies.
}

}
