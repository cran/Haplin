\name{hapSim}

\alias{hapSim}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Simulation of genetic data in Haplin format}

\description{Simulates genetic data in Haplin format, consisting of fetal effects, maternal effects and/or parent-of-origin effects. 
Allows for simulation of both autosomal and X-linked markers, assuming Hardy-Weinberg equilibrium.
Enables stratified simulations for gene-environment interaction analyses, i.e the input (relative risks, number of cases etc) may vary across different strata.}

\usage{
hapSim(nall, n.strata = 1, cases, controls, haplo.freq,
RR, RRcm, RRcf, RRstar, RR.mat, RRstar.mat, 
gen.missing.cases = NULL, gen.missing.controls = NULL,
n.sim = 1000, xchrom = F, sim.comb.sex = "double", BR.girls, dire = "simfiles", 
ask = TRUE, verbose = TRUE, cpus = 1)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
	\item{nall}{A vector of the number of alleles at each locus.}
	\item{n.strata}{The number of strata.}
	\item{cases}{A list of the number of case families. Each element is a vector of the number of families of the specified family design(s) in the corresponding stratum. The possible family designs, i.e., the possible names of the elements, are \code{"mfc"} (full triad), \code{"mc"} (mother-child dyad), \code{"fc"} (father-child dyad) or \code{"c"} (a single case child). See Details for a thorough description.}
	\item{controls}{A list of the number of control families. Each element is a vector of the number of families of the specified family design(s) in the corresponding stratum. The possible family designs are \code{"mfc"} (full triad), \code{"mc"} (mother-child-dyad), \code{"fc"} (father-child dyad), \code{"mf"} (mother-father dyad), \code{"c"} (a single control child), \code{"m"} (a single control mother) or \code{"f"} (a single control father). See Details for a thorough description.}
	\item{haplo.freq}{A list of which each element is a numeric vector of the haplotype frequencies in each stratum. The frequencies will be normalized so that they sum to one. The Details section shows how to implement this argument in agreement with the possible haplotypes.}
	\item{RR}{A list of which each element is a numeric vector of the relative risks in each stratum. The Details section shows how to implement this argument in agreement with the possible haplotypes.}
	\item{RRcm}{A list of numeric vectors. Each vector contains the relative risks associated with the haplotypes transmitted from the mother for this stratum. See Details for description of how to implement this argument in agreement with the possible haplotypes.}
	\item{RRcf}{A list of numeric vectors. Each vector contains the relative risks associated with the haplotypes transmitted from the father for this stratum. See Details for description of how to implement this argument in agreement with the possible haplotypes.}
	\item{RRstar}{A list of numeric vectors. Estimates how much double-dose children would deviate from the risk expected in a multiplicative dose-response relationship.}
	\item{RR.mat}{The interpretation is similar to \code{RR} when simulating genetic data with maternal effects.}
	\item{RRstar.mat}{The interpretation is parallel to \code{RRstar} when simulating genetic data with maternal effects.}
	\item{gen.missing.cases}{Generates missing values at random for the case families. Set to \code{NULL} by default, i.e., no missing values generated. See Details for description of how to implement this argument.}
	\item{gen.missing.controls}{Generates missing values at random for the control families. Set to \code{NULL} by default, i.e., no missing values generated. See Details for description of how to implement this argument.}
	\item{n.sim}{The number of simulations, i.e., the number of simulated data files.}
	\item{xchrom}{Logical. Equals \code{FALSE} by default, which indicates simulation of autosomal markers. If \code{TRUE}, \code{hapSim} simulates X-linked genes.}
	\item{sim.comb.sex}{To be used with \code{xchrom = TRUE}. A character value which specifies how to handle gender differences on the X-chromosome. If \code{"single"}, the effect of a (single) allele in males is equal to the effect of a single allele dose in females, and similarly, if \code{"double"}, a single allele in males has the same effect as a double allele dose in females. Default is \code{"double"}, which corresponds to X-inactivation.} 
	\item{BR.girls}{To be used with \code{xchrom = TRUE}. Gives the ratio of baseline risk for females to the baseline risk for males.}
	\item{dire}{Gives the directory of the simulated data files.}
	\item{ask}{Logical. If \code{TRUE}, \code{hapSim} will ask before overwriting the files in an already existing directory.}
	\item{verbose}{Logical. Default is \code{TRUE}, which means that the file name is displayed for each iteration. Works only when \code{cpus = 1}.}
	\item{cpus}{Allows simulations to be performed in parallel. The \code{cpus} argument should preferably be set to the number of available cores. If set lower, it will save some capacity for other processes to run. Setting it too high should not cause any serious problems.}
}

\details{\code{hapSim} simulates allele values for case and control families at multiple markers (typically in LD) simultaneously. The number of markers/SNPs involved will typically be in the range 1 to 6. Data are simulated to produce relative risks of disease as specified by the user input. Simulations can be performed separately over a number of strata so as to simulate gene-environment interactions.

Specifying haplotype risks:\cr
The number of haplotypes used in the simulations is determined by the \code{nall} argument, since \code{prod(nall)} different haplotypes can be made from the specified number of markers, \code{length(nall)}. The arguments \code{haplo.freq}, \code{RR}, \code{RRcm}, \code{RRcf}, \code{RRstar}, \code{RR.mat}, and \code{RRstar.mat} are all lists where each element represents a stratum. Within each stratum, the arguments are vectors of length equal to the number of haplotypes, specifying the relative risk etc. associated with each haplotype.
The stratum specific arguments may be simplified if the number of strata is one, or if the arguments are equal across all strata.

The haplotypes are determined by creating all possible haplotypes from the given markers, in a sequence where the first marker varies mostly quickly. For instance, if \code{nall = c(3,2)}, the first marker has 3 alleles, the second has 2, and 6 haplotypes are possible. Taken in order, the haplotypes are 1-1, 2-1, 3-1, 1-2, 2-2, and 3-2. When specifying, say, \code{RR = c(1,2,1,1,1,1)} the haplotype 2-1 has a double risk compared to the rest. With, for instance, two strata, the specification \code{RR = list(c(1,2,1,1,1,1), c(1,1,1,1,1,1))} would mean that the risk associated with 2-1 is elevated only in the first stratum, not the second. 

The simplest example would be with \code{nall = c(2)} and \code{RR = c(1,2)}, which would simulate a single SNP where the second allele has a double risk.

Output file format:\cr
The format of the simulated files is relatively flexible and allows multi-allelic markers and various designs.
If both case and control families are present, the simulated files contain a leading column of the case/control status (1/0). 
If \code{xchrom=TRUE}, the neighboring column to the left of the genetic data contains the sex information (1 = male, 2 = female). 
Each line represents genotypes for a case or control triad. 

There are six columns for each locus, two for the mother (M), two for the father (F) and two for the child (C). The columns are placed in the following sequence: M11 M12 F11 F12 C11 C12 M21 M22 F21 F22 C21 C22... etc, 
where the first number indicates marker, and the second number indicates the first or second allele at this locus. Columns are separated by a single white space, and missing data are coded as NA.

Intermediate designs, for instance mother-child dyads, are represented as full triads with columns of absent family members set to missing. In the case of a pure case-control design, however, each line represents a single individual, and there are no columns representing mothers and fathers.\cr

There are no row or column names in the files.

Some examples are given below. See \url{https://people.uib.no/gjessing/genetics/software/haplin/docu/haplin_data_format.html} for a thorough description of the Haplin format. Note that this description separates the two alleles for an individual within a locus by a semi-colon, such as 1;2. This is, however, not necessary.\cr
Confer the document \url{https://people.uib.no/gjessing/genetics/software/haplin/docu/Haplin_power.pdf} for details and examples on how to perform the simulations.

\code{gen.missing.cases} and \code{gen.missing.controls} are flexible arguments. By default, both equal NULL, which means that no missing data are generated at random.
If the arguments are single numbers, missing data are generated at random with this proportion for all cases and/or controls. 
If the arguments are vectors of length equal to the number of loci, missing data are generated with the corresponding proportion for each locus.
The arguments can also be matrices with the number of rows equal to the number of loci and three columns. 
Each row corresponds to a locus, and the columns correspond to mothers, fathers and children, respectively.}


\references{
Web Site: \url{https://people.uib.no/gjessing/genetics/software/haplin/}
}

\author{Miriam Gjerdevik,\cr
with Hakon K. Gjessing\cr
Professor of Biostatistics\cr
Division of Epidemiology\cr
Norwegian Institute of Public Health\cr

\email{hakon.gjessing@uib.no}}

\seealso{\code{\link[Haplin]{haplin}}, \code{\link[Haplin]{hapRun}}, \code{\link[Haplin]{hapPower}}}


\examples{
\dontrun{
## Simulate genetic data (100 files) at two diallelic markers, consisting of fetal effects 
## corresponding to haplo.freq = rep(0.25, 4), RR = c(2,1,1,1) and RRstar = c(1,1,1,1),
## for the combination of 1000 case triads and 1000 control triads with no missing data.
## Only one stratum.
hapSim(nall = c(2,2), n.strata = 1, cases = c(mfc=1000),
controls = c(mfc=1000), haplo.freq = rep(0.25,4),
RR = c(2,1,1,1), RRstar = c(1,1,1,1), n.sim = 100, dire = "simfiles")

## Simulate genetic data (100 files) at two diallelic markers,
## consisting of fetal and maternal effects corresponding to 
## haplo.freq = rep(0.25, 4), RR = c(2,1,1,1), RRstar = c(1,1,1,1),
## RR.mat = c(2,1,1,1) and RRstar.mat = c(1,1,1,1),
## for 1000 case triads and zero control families. 
## One percent of the case triads are missing at random. One stratum only.
hapSim(nall = c(2,2), n.strata=1, cases = c(mfc=1000),
controls = c(mfc=0), haplo.freq = rep(0.25,4), RR = c(2,1,1,1),
RRstar = c(1,1,1,1), RR.mat = c(2,1,1,1), RRstar.mat = c(1,1,1,1), 
gen.missing.cases = 0.01, n.sim = 100, dire = "simfiles")

## Simulate genetic data (100 files) at two diallelic markers. In the first stratum, 
## we have a combination of 500 case triads and 500 control triads with
## haplo.freq = rep(0.25, 4), RR = c(2,1,1,1) and RRstar = c(1,1,1,1).
## In the second stratum, we have 300 case triads and 500 control triads with
## haplo.freq = rep(0.25, 4), RR = c(1,1,1,1) and RRstar = c(1,1,1,1).
## One percent of the control triads are missing at random in the first stratum.
hapSim(nall = c(2,2), n.strata= 2, cases = list(c(mfc=500),c(mfc=300)), 
controls = c(mfc=500),haplo.freq = rep(0.25,4),
RR = list(c(2,1,1,1),c(1,1,1,1)), RRstar = c(1,1,1,1),
gen.missing.controls = list(0.01,NULL), n.sim = 100, dire = "simfiles")
}
}
