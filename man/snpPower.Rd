\name{snpPower}

\alias{snpPower}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Power calculations for a single SNP}

\description{Calculates power for a single SNP. Allows for power computations of several combinations simultaneously.}

\usage{
snpPower(cases, controls, RR, MAF, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
	\item{cases}{A list of the number of case families. Each element contains the number of families of a specified family design. The possible family designs, i.e. the possible names of the elements, are \code{"mfc"} (full triad), \code{"mc"} (mother-child dyad), \code{"fc"} (father-child dyad) or \code{"c"} (a single case child).}
	\item{controls}{A list of the number of control families. Each element contains the number of families of a specified family design. The possible family designs are \code{"mfc"} (full triad), \code{"mc"} (mother-child-dyad), \code{"fc"} (father-child dyad), \code{"mf"} (mother-father dyad), \code{"c"} (a single control child), \code{"m"} (a single control mother) or \code{"f"} (a single control father).}
	\item{RR}{A numeric vector of the relative risks (the effect sizes of interest).}
	\item{MAF}{A numeric vector of the minor allele frequencies.}
	\item{alpha}{A numeric vector of the Type I Errors. Equals 0.05 by default.}
}

\details{\code{snpPower} computes power for a single SNP by counting the number of "real" case alleles, "real" control alleles and pseudo-control alleles. 
The pseudo-control alleles are the non-transmitted alleles, possibly from both case families and control families. It assumes a multiplicative dose-response model.
\code{snpPower} uses the asymptotic normal approximation for the natural logarithm of the odds ratio for calculating power 
(the relative risks and odds ratios are used interchangeably due to the "rare disease assumption").\cr

\code{snpPower} allows for power calculations for mixtures of the possible case family designs and control family designs. 
The argument \code{cases} could, for example, consist of a combination of 1000 full case triads (family design \code{"mfc"}) and 500 single case children (family design \code{"c"}). 
It is also feasible to compute power for several combinations of the input variables simultaneously. See Examples for further details.
}

\value{\code{snpPower} returns a data frame containing the combinations of input variables and the corresponding power calculations.}

\references{
Skare O, Jugessur A, Lie RT, Wilcox AJ, Murray JC, Lunde A, Nguyen TT, Gjessing HK. Application of a novel hybrid study design to explore gene-environment interactions in orofacial clefts. Annals of Human Genetics (2012) 76, pp. 221-236.\cr

Web Site: \url{http://folk.uib.no/gjessing/genetics/software/haplin/}
}

\author{Miriam Gjerdevik,\cr
with Hakon K. Gjessing\cr
Professor of Biostatistics\cr
Division of Epidemiology\cr
Norwegian Institute of Public Health\cr

\email{hakon.gjessing@uib.no}}

\seealso{\code{\link[Haplin]{haplin}}, \code{\link[Haplin]{snpSampleSize}}}

\examples{
\dontrun{
## Compute power for a single SNP,
## for the combination of 1000 case triads, RR = 1.2, MAF = 0.1 and alpha = 0.05
snpPower(cases = list(mfc = 1000), controls = list(mfc = 0), RR = 1.2, MAF = 0.1)

## Compute power for a single SNP,
## for the combination of 1000 case triads and 500 single case children (altogether),
## 5000 control triads, RR = 1.1, MAF = 0.1 and alpha = 0.05
snpPower(cases = list(mfc = 1000, c = 500), controls = list(mfc=5000),
RR = 1.1, MAF = 0.1, alpha = 0.05)

## Compute power for a single SNP,
## for the combination of 500 case triads, 10000 control triads, 
## relative risk of 1.2 and minor allele frequency of 0.1,
## and also for the combination of 1000 case triads, 10000 control triads,
## relative risk of 1.1 and minor allele frequency of 0.1
snpPower(cases = list(mfc = c(500, 1000)), controls = list(mfc = 10000),
RR = c(1.2, 1.1), MAF = 0.1)
}
}

