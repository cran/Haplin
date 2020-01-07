\name{hapPowerAsymp}

\alias{hapPowerAsymp}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Asymptotic power calculations for genetic association analyses with Haplin}

\description{Computes the asymptotic power for genetic analyses assessing fetal effects, maternal effects and/or parent-of-origin effects. Effects of X-chromosome genes and gene-environment interaction effects are also allowed. 
}

\usage{
hapPowerAsymp(nall = 2, n.strata = 1, cases, controls, haplo.freq,
RR, RRcm, RRcf, RRstar, RR.mat, RRstar.mat,
xchrom = F, sim.comb.sex = "double", BR.girls, 
response = "mult", alpha = 0.05, ...)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
	\item{nall}{A vector of the number of alleles at each locus. By default a diallellic SNP.}
	\item{n.strata}{The number of strata.}
	\item{cases}{A list of the number of case families. Each element is a vector of the number of families of the specified family design in the corresponding stratum. The possible family designs, i.e., the possible names of the elements, are \code{"mfc"} (full triad), \code{"mc"} (mother-child dyad), \code{"fc"} (father-child dyad) or \code{"c"} (a single case child).}
	\item{controls}{A list of the number of control families. Each element is a vector of the number of families of the specified family design in the corresponding stratum. The possible family designs are \code{"mfc"} (full triad), \code{"mc"} (mother-child-dyad), \code{"fc"} (father-child dyad), \code{"mf"} (mother-father dyad), \code{"c"} (a single control child), \code{"m"} (a single control mother) or \code{"f"} (a single control father).}
	\item{haplo.freq}{A list of which each element is a numeric vector of the haplotype frequencies in each stratum. The frequencies are normalized and sum to one. The Details section shows how to implement this argument in agreement with the possible haplotypes.}
	\item{RR}{A list of which each element is a numeric vector of the relative risks in each stratum. The Details section shows how to implement this argument in agreement with the possible haplotypes.}
	\item{RRcm}{A list of numeric vectors. Each vector contains the relative risks associated with the haplotypes transmitted from the mother for this stratum. See Details for description of how to implement this argument in agreement with the possible haplotypes.}
	\item{RRcf}{A list of numeric vectors. Each vector contains the relative risks associated with the haplotypes transmitted from the father for this stratum. See Details for description of how to implement this argument in agreement with the possible haplotypes.}
	\item{RRstar}{A list of numeric vectors. Estimates how much double-dose children would deviate from the risk expected in a multiplicative dose-response relationship.}
	\item{RR.mat}{The interpretation is similar to \code{RR} but for maternal genetic effects.}
	\item{RRstar.mat}{The interpretation is similar to \code{RRstar} but for maternal genetic effects.}
	\item{xchrom}{Logical. Equals \code{FALSE} by default, which indicates analyses of autosomal markers. If \code{TRUE}, \code{hapPowerAsymp} performs power analyses of X-linked markers.}
	\item{sim.comb.sex}{To be used with \code{xchrom = TRUE}. A character value that specifies how to handle gender differences on the X-chromosome. If \code{"single"}, the effect of a (single) allele in males is equal to the effect of a single allele dose in females, and similarly, if \code{"double"}, a single allele in males has the same effect as a double allele dose in females. Default is "double", which corresponds to X-inactivation.} 
	\item{BR.girls}{To be used with \code{xchrom = TRUE}. Gives the ratio of baseline risk for females relative to the baseline risk for males.}
	\item{response}{The default value "mult" specifies a multiplicative dose-response model. \code{response = "free"} is not yet implemented.}
	\item{alpha}{\code{alpha} is the Type I Error probablity. Equals 0.05 by default.}
	\item{...}{Could include argument \code{reference}. By default, the most frequent allele or haplotype is chosen as reference. The reference haplotype can be set explicitly by giving a numeric value for the reference argument.} 
	}


\details{The Haplin framework includes different modules for assessing genetic effects: \code{\link[Haplin]{haplin}}, 
\code{\link[Haplin]{haplinStrat}} and \code{\link[Haplin]{haplinSlide}}. \code{hapPowerAsymp} computes the asymptotic power for these analyses, which enables power calculations of fetal effects, maternal effects and/or parent-of-origin effects. Various family designs, i.e., triads, case-control, the hybrid design, and all intermediate designs, are possible.
It also allows power calculation of gene-environment interaction effects and effects on X-chromosome markers.

\code{\link[Haplin]{hapPower}} computes power through "brute force" simulations using \code{\link[Haplin]{hapRun}}. This is a robust way of checking software implementations, asymptotic approximations
and attained significance level. However, both power and the corresponding sample size calculations can be performed much more efficiently using asymptotic approximations.
The asymptotic power is calculated applying the non-centrality parameter of the Wald tests, which use the asymptotic normal distribution of the log-scale parameter. The function \code{hapCovar} (used by \code{hapPowerAsymp}) computes the variance-covariance matrix 
by applying the log-linear model combined with the EM algorithm.

Specifying haplotype risks:\cr
The number of haplotypes used in the simulations is determined by the \code{nall} argument, since \code{prod(nall)} different haplotypes can be made from the specified number of markers, \code{length(nall)}. The arguments \code{haplo.freq}, \code{RR}, \code{RRcm}, \code{RRcf}, \code{RRstar}, \code{RR.mat}, and \code{RRstar.mat} are all lists where each element represents a stratum. Within each stratum, the arguments are vectors of length equal to the number of haplotypes, specifying the relative risk etc. associated with each haplotype.
The stratum specific arguments may be simplified if the number of strata is one, or if the arguments are equal across all strata.
The haplotypes are determined by creating all possible haplotypes from the given markers, in a sequence where the first marker varies mostly quickly. For instance, if \code{nall = c(3,2)}, the first marker has 3 alleles, the second has 2, and 6 haplotypes are possible. Taken in order, the haplotypes are 1-1, 2-1, 3-1, 1-2, 2-2, and 3-2. When specifying, say, \code{RR = c(1,2,1,1,1,1)} the haplotype 2-1 has a double risk compared to the rest. With, for instance, two strata, the specification \code{RR = list(c(1,2,1,1,1,1), c(1,1,1,1,1,1))} would mean that the risk associated with 2-1 is elevated only in the first stratum, not the second.
The simplest example would be with \code{nall = c(2)} and \code{RR = c(1,2)}, which would simulate a single SNP where the second allele has a double risk.

Specifying genetic effects:\cr
Standard fetal effects are specified by the arguments \code{RR} and \code{RRstar}, whereas parent-of-origin effects are addressed by the arguments \code{RRcm}, \code{RRcf} and \code{RRstar}.  Maternal effects are included by the additional arguments \code{RRmat} and \code{RRstar.mat}.

}

\value{\code{hapPowerAsymp} returns the asymptotic power for the relevant genetic effects. The first element of the list depicts the power for each haplotype analyzed separately. 
If there are more than two possible haplotypes, the second element displays the overall power for all haplotypes combined.}

\references{
Web Site: \url{http://folk.uib.no/gjessing/genetics/software/haplin/}
}

\author{Miriam Gjerdevik,\cr
with Hakon K. Gjessing\cr
Professor of Biostatistics\cr
Division of Epidemiology\cr
Norwegian Institute of Public Health\cr

\email{hakon.gjessing@uib.no}}

\seealso{\code{\link[Haplin]{haplin}}, \code{\link[Haplin]{haplinSlide}},  \code{\link[Haplin]{haplinStrat}}, \code{\link[Haplin]{hapSim}}, \code{\link[Haplin]{hapRun}}, \code{\link[Haplin]{snpPower}}, \code{\link[Haplin]{hapPower}}}


\examples{

## Calculate the asymptotic power for a triad design 
## when the minor allele increases the fetal risk by twofold.
## Assumes a multiplicative dose-response relationship.
hapPowerAsymp(nall = c(2), n.strata = 1, cases = list(c(mfc=120)), 
haplo.freq = c(0.1,0.9), RR = c(2,1), RRstar = c(1,1))

## Calculate the asymptotic power for the hybrid design when 
## the minor allele increases the fetal risk by twofold 
## in the first stratum and no effect is seen in the second
## i.e., gene-environment interaction (GxE) effects.
hapPowerAsymp(nall = c(2), n.strata = 2, cases = list(c(mfc=100)), 
controls = list(c(mfc=100)), haplo.freq = c(0.1,0.9), 
RR = list(c(2,1), c(1,1)), RRstar = c(1,1))

## Calculate the asymptotic GxE power assessing maternal- and
## parent-of-origin effects at two diallelic loci.
hapPowerAsymp(nall = c(2,2), n.strata = 2, cases = c(mfc=500), 
haplo.freq = c(0.1,0.2,0.3,0.4), 
RRcm = list(c(3,1,1,1), c(1,1,1,1)), RRcf = c(1,1,1,1), RRstar = c(1,1,1,1),
RR.mat = list(c(1.5,1,1,1),c(1,1,1,1)), RRstar.mat = c(1,1,1,1))

}

