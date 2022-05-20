\name{snpSampleSize}

\alias{snpSampleSize}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Sample size calculations for a single SNP}

\description{Sample size calculations for a single SNP. Allows for sample size calculations of several combinations simultaneously.}

\usage{
snpSampleSize(fam.cases, fam.controls, fraction = 0.5,
RR, MAF, alpha = 0.05, power = 0.80)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
	\item{fam.cases}{A character vector of the case family design. The possible family designs are \code{"mfc"} (full triad), \code{"mc"} (mother-child dyad), \code{"fc"} (father-child dyad) and \code{"c"} (a single case child).}
	\item{fam.controls}{A character vector of the control family design. The possible family designs are \code{"mfc"} (full triad), \code{"mc"} (mother-child-dyad), \code{"fc"} (father-child dyad), \code{"mf"} (mother-father dyad), \code{"c"} (a single control child), \code{"m"} (a single control mother), \code{"f"} (a single control father) or \code{"no_controls"} (no control families).}
	\item{fraction}{A numeric vector of the proportion of case families. Equals 0.5 by default, i.e. there are as many case families as control families. If \code{fam.controls} equals \code{"no_controls"}, \code{fraction} is automatically set to 1.}
	\item{RR}{A numeric vector of the relative risks (the effect sizes of interest).}
	\item{MAF}{A numeric vector of the minor allele frequencies.}
	\item{alpha}{A numeric vector of the Type I Errors. Equals 0.05 by default.}
	\item{power}{A numeric vector of the desired probability of identifying a difference in the relative risks. Default is 0.80.}
}

\details{\code{snpSampleSize} computes the number of case and control families required for a single SNP to attain the desired power. 
It assumes a multiplicative dose-response model.
\code{snpSampleSize} calculates the fraction of case alleles corresponding to the given family designs 
and then uses the asymptotic normal approximation for the natural logarithm of the odds ratio for calculating the sample sizes 
(the relative risks and odds ratios are used interchangeably due to the "rare disease assumption").\cr

\code{snpSampleSize} allows for sample size calculations of several combinations of the input variables at once. The Examples section provides further details. 
}

\value{\code{snpSampleSize} returns a data frame containing the combinations of input variables and the corresponding sample size calculations.}

\references{
Skare O, Jugessur A, Lie RT, Wilcox AJ, Murray JC, Lunde A, Nguyen TT, Gjessing HK. Application of a novel hybrid study design to explore gene-environment interactions in orofacial clefts. Annals of Human Genetics (2012) 76, pp. 221-236.\cr

Web Site: \url{https://people.uib.no/gjessing/genetics/software/haplin/}
}

\author{Miriam Gjerdevik,\cr
with Hakon K. Gjessing\cr
Professor of Biostatistics\cr
Division of Epidemiology\cr
Norwegian Institute of Public Health\cr

\email{hakon.gjessing@uib.no}}

\seealso{\code{\link[Haplin]{haplin}}, \code{\link[Haplin]{snpPower}}}

\examples{

## Compute sample sizes for a single SNP,
## when the specified case family design is a full triad, there are no control families,
## RR = 1.1, MAF = 0.1, alpha = 0.05 and power = 0.9
snpSampleSize(fam.cases = "mfc", fam.controls = "no_controls", 
RR = 1.1, MAF = 0.1, alpha = 0.05, power = 0.9)

## Compute sample sizes for a single SNP,
## for the combination of case triads, control triads, fraction = 0.5, RR = 1.2, 
## MAF = 0.1, alpha = 0.05 and power = 0.8, and also for the combination of case triads,
## control children, fraction = 0.5, RR = 1.2, MAF = 0.1, alpha = 0.05 and power = 0.9
snpSampleSize(fam.cases = "mfc", fam.controls = c("mfc", "c"), 
RR = 1.2, MAF = 0.1, alpha = 0.05, power = c(0.8, 0.9))

}

