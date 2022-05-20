\name{hapRelEff}

\alias{hapRelEff}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Relative efficiency comparing different study designs in genetic association analysis with Haplin}

\description{Computes the relative efficiency for different study designs in genetic association analysis. 
}

\usage{
hapRelEff(nall = 2, cases.comp, controls.comp, 
cases.ref, controls.ref, haplo.freq,
RR, RRcm, RRcf, RRstar, RR.mat, RRstar.mat,
xchrom = F, sim.comb.sex = "double", BR.girls, 
response = "mult", ...)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
	\item{nall}{A vector of the number of alleles at each locus. By default a diallellic SNP.}
	\item{cases.comp}{A list of the number of case families in the comparison design. Its element is a vector of the number of families of the specified family design. The possible family designs, i.e., the possible names of the elements, are \code{"mfc"} (full triad), \code{"mc"} (mother-child dyad), \code{"fc"} (father-child dyad) or \code{"c"} (a single case child).}		
	\item{controls.comp}{A list of the number of control families in the comparison design. Its element is a vector of the number of families of the specified family design. The possible family designs are \code{"mfc"} (full triad), \code{"mc"} (mother-child-dyad), \code{"fc"} (father-child dyad), \code{"mf"} (mother-father dyad), \code{"c"} (a single control child), \code{"m"} (a single control mother) or \code{"f"} (a single control father).}
	\item{cases.ref}{A list of the number of case families in the reference design. The options are the same as for \code{cases.comp}.}
	\item{controls.ref}{A list of the number of control families in the reference design. The options are the same as for \code{controls.comp}.}
	\item{haplo.freq}{A list of which its element is a numeric vector of the haplotype frequencies. The frequencies are normalized and sum to one. The Details section shows how to implement this argument in agreement with the possible haplotypes.}
	\item{RR}{A list of which its element is a numeric vector of the relative risks. The Details section shows how to implement this argument in agreement with the possible haplotypes.}
	\item{RRcm}{A numeric vector in list format containing the relative risks associated with the haplotypes transmitted from the mother. See Details for description of how to implement this argument in agreement with the possible haplotypes.}
	\item{RRcf}{A numeric vector in list format containing the relative risks associated with the haplotypes transmitted from the father. See Details for description of how to implement this argument in agreement with the possible haplotypes.}
	\item{RRstar}{A numeric vector in list format. Estimates how much double-dose children would deviate from the risk expected in a multiplicative dose-response relationship.}
	\item{RR.mat}{The interpretation is similar to \code{RR} but for maternal genetic effects.}
	\item{RRstar.mat}{The interpretation is similar to \code{RRstar} but for maternal genetic effects.}
	\item{xchrom}{Logical. Equals \code{FALSE} by default, which indicates analyses of autosomal markers. If \code{TRUE}, analyses are perfomed on X-linked markers.}
	\item{sim.comb.sex}{To be used with \code{xchrom = TRUE}. A character value that specifies how to handle gender differences on the X-chromosome. If \code{"single"}, the effect of a (single) allele in males is equal to the effect of a single allele dose in females, and similarly, if \code{"double"}, a single allele in males has the same effect as a double allele dose in females. Default is "double", which corresponds to X-inactivation.} 
	\item{BR.girls}{To be used with \code{xchrom = TRUE}. Gives the ratio of baseline risk for females relative to the baseline risk for males.}
	\item{response}{The default value "mult" specifies a multiplicative dose-response model. \code{response = "free"} is not yet implemented.}
	\item{...}{Could include argument \code{reference}. By default, the most frequent allele or haplotype is chosen as reference. The reference haplotype can be set explicitly by giving a numeric value for the reference argument.} 
	}


\details{\code{hapRelEff} compares two study designs for genetic association analysis, using a term called relative efficiency. 
The relative efficiency is defined as the ratio of variances of estimators for the same parameter, computed from two
different designs, or equivalently, the ratio of the sample sizes needed for each of the two designs to achieve the same level and power. The number of genotyped individuals within each design is accounted for.
The relative efficiency estimated under the null hypothesis, i.e., when all relative risks are equal to one, is known as the Pitman efficiency.
To compute the variance-covariance estimates for each design, \code{hapRelEff} calls the function \code{hapCovar}, which calculates the asymptotic variance-covariance matrix 
by applying the log-linear model combined with the EM algorithm.

The relative efficiency can be computed for fetal effects, maternal effects and parent-of-origin effect, 
as well as effects on X-chromosome markers. 
Various study designs, i.e., case-parent triads, the standard case-control design, the hybrid design, and all intermediate designs, can be compared.

Note that the exact number of case families and control families in the reference or comparison design is irrelevant, as this will be accounted for in the relative efficiency estimate.
However, the ratio of control families to case families within the reference or comparison design must be specified correctly. See the Examples section.

Specifying haplotype risks:\cr
The number of haplotypes used in the simulations is determined by the \code{nall} argument, since \code{prod(nall)} different haplotypes can be made from the specified number of markers, \code{length(nall)}. 
The arguments \code{haplo.freq}, \code{RR}, \code{RRcm}, \code{RRcf}, \code{RRstar}, \code{RR.mat}, and \code{RRstar.mat} are all lists, containing vectors of length equal to the number of haplotypes, specifying the relative risk etc. associated with each haplotype.
However, the function will work without using the list format.
The haplotypes are determined by creating all possible haplotypes from the given markers, in a sequence where the first marker varies mostly quickly. 
For instance, if \code{nall = c(3,2)}, the first marker has 3 alleles, the second has 2, and 6 haplotypes are possible. Taken in order, the haplotypes are 1-1, 2-1, 3-1, 1-2, 2-2, and 3-2. 
When specifying, say, \code{RR = c(1,2,1,1,1,1)} the haplotype 2-1 has a double risk compared to the rest. 
The simplest example would be with \code{nall = c(2)} and \code{RR = c(1,2)}, which would simulate a single SNP where the second allele has a double risk.

Specifying genetic effects:\cr
Standard fetal effects are specified by the arguments \code{RR} and \code{RRstar}, 
whereas parent-of-origin effects are addressed by the arguments \code{RRcm}, \code{RRcf} and \code{RRstar}. 
Maternal effects are included by the additional arguments \code{RRmat} and \code{RRstar.mat}.

}

\value{\code{hapRelEff} returns the relative efficiency estimate, comparing two study designs. The first element of the list depicts the relative efficiency for each haplotype analyzed separately. 
If there are more than two possible haplotypes, the second element displays the overall relative efficiency for all haplotypes combined.}

\references{
Web Site: \url{https://people.uib.no/gjessing/genetics/software/haplin/}
}

\author{Miriam Gjerdevik,\cr
with Hakon K. Gjessing\cr
Professor of Biostatistics\cr
Division of Epidemiology\cr
Norwegian Institute of Public Health\cr

\email{hakon.gjessing@uib.no}}

\seealso{\code{\link[Haplin]{haplin}}, \code{\link[Haplin]{hapPowerAsymp}}, \code{\link[Haplin]{hapRun}}, \code{\link[Haplin]{snpPower}}, \code{\link[Haplin]{hapPower}}}


\examples{

## Child effects: Calculate the efficiency of the standard case-control design 
## (with an equal number of case and control children) 
## relative to the case-parent triad design
## under the null hypothesis when the minor allele frequency is 0.1.
hapRelEff(nall = c(2), cases.comp = c(c=1), 
controls.comp = c(c=1), cases.ref = c(mfc=1),
haplo.freq = c(0.1,0.9), RR = c(1,1))
	
## Child effects: Calculate the efficiency of the standard case-control design, 
## with twice as many cases as controls, relative to the case-parent triad design
## under the null hypothesis when the minor allele frequency is 0.2.
hapRelEff(nall = c(2), cases.comp = c(c=2), 
controls.comp = c(c=1), cases.ref = c(mfc=1),
haplo.freq = c(0.2,0.8), RR = c(1,1))

## Child and maternal effects: Calculate the efficiency of the case-mother dyad design 
## relative to the case-parent triad design
## under the null hypothesis when the minor allele frequency is 0.1.
hapRelEff(nall = c(2), cases.comp = c(mc=1), cases.ref = c(mfc=1),
haplo.freq = c(0.1,0.9), RR = c(1,1), RR.mat=c(1,1))

## PoO effects: Calculate the efficiency of the full hybrid design, 
## with twice as many control families as case families, 
## relative to the case-parent triad design
## under the null hypothesis when the minor allele frequency is 0.1.
hapRelEff(nall = c(2), cases.comp = c(mfc=1), 
controls.comp = c(mfc=2), cases.ref = c(mfc=1),
haplo.freq = c(0.1,0.9), RRcm = c(1,1), RRcf = c(1,1))

}

