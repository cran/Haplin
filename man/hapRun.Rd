\name{hapRun}

\alias{hapRun}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Simulates genetic data and runs Haplin for each simulation}

\description{Calculates Haplin results by first simulating genetic data, allowing a various number of family designs, and then running Haplin on the simulations.
The simulated data may contain of fetal effects, maternal effects and/or parent-of-origin effects.
The function allows for simulations and calculations on both autosomal and X-chromosome markers, 
assuming Hardy-Weinberg equilibrium.
It enables simulation and calculation of gene-environment interaction effects, i.e, the input (relative risks, number of cases etc.) may vary across strata.
\code{hapRun} calls \code{\link[Haplin]{haplin}}, \code{\link[Haplin]{haplinStrat}} or \code{\link[Haplin]{haplinSlide}} to run on the simulated data files.
}

\usage{
hapRun(nall, n.strata= 1, cases, controls, haplo.freq,
RR, RRcm, RRcf, RRstar, RR.mat, RRstar.mat, hapfunc = "haplin",
gen.missing.cases = NULL, gen.missing.controls = NULL, 
n.sim = 1000, xchrom = FALSE, sim.comb.sex = "double", BR.girls, dire,
ask = TRUE, cpus = 1, slaveOutfile = "", ...)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
	\item{nall}{A vector of the number of alleles at each locus.}
	\item{n.strata}{The number of strata.}
	\item{cases}{A list of the number of case families. Each element is a vector of the number of families of the specified family design(s) in the corresponding stratum. The possible family designs, i.e., the possible names of the elements, are \code{"mfc"} (full triad), \code{"mc"} (mother-child dyad), \code{"fc"} (father-child dyad) or \code{"c"} (a single case child). See Details for a thorough description.}
	\item{controls}{A list of the number of control families. Each element is a vector of the number of families of the specified family design(s) in the corresponding stratum. The possible family designs are \code{"mfc"} (full triad), \code{"mc"} (mother-child-dyad), \code{"fc"} (father-child dyad), \code{"mf"} (mother-father dyad), \code{"c"} (a single control child), \code{"m"} (a single control mother) or \code{"f"} (a single control father). See Details for a thorough description.}
	\item{haplo.freq}{A list of which each element is a numeric vector of the haplotype frequencies in each stratum. The frequencies are normalized and sum to one. The Details section shows how to implement this argument in agreement with the possible haplotypes.}
	\item{RR}{A list of which each element is a numeric vector of the relative risks in each stratum. The Details section shows how to implement this argument in agreement with the possible haplotypes.}
	\item{RRcm}{A list of numeric vectors. Each vector contains the relative risks associated with the haplotypes transmitted from the mother for this stratum. See Details for description of how to implement this argument in agreement with the possible haplotypes.}
	\item{RRcf}{A list of numeric vectors. Each vector contains the relative risks associated with the haplotypes transmitted from the father for this stratum. See Details for description of how to implement this argument in agreement with the possible haplotypes.}
	\item{RRstar}{A list of numeric vectors. Estimates how much double-dose children would deviate from the risk expected in a multiplicative dose-response relationship.}
	\item{RR.mat}{The interpretation is similar to \code{RR} when simulating genetic data with maternal effects.}
	\item{RRstar.mat}{The interpretation is similar to \code{RRstar} when simulating genetic data with maternal effects.}
	\item{hapfunc}{Defines which haplin function to run, the options being \code{"haplin"}, \code{"haplinSlide"} or \code{"haplinStrat"}. \code{"haplinSlide"} is however only partially implemented.}
	\item{gen.missing.cases}{Generates missing values at random for the case families. Set to \code{NULL} by default, i.e., no missing values generated. See Details for description of how to implement this argument.}
	\item{gen.missing.controls}{Generates missing values at random for the control families. Set to \code{NULL} by default, i.e., no missing values generated. See Details for description of how to implement this argument.}
	\item{n.sim}{The number of simulations, i.e., the number of simulated data files.}
	\item{xchrom}{Logical. Equals \code{FALSE} by default, which indicates simulation of autosomal markers. If \code{TRUE}, hapSim simulates X-linked genes.}
	\item{sim.comb.sex}{To be used with \code{xchrom = TRUE}. A character value that specifies how to handle gender differences on the X-chromosome. If \code{"single"}, the effect of a (single) allele in males is equal to the effect of a single allele dose in females, and similarly, if \code{"double"}, a single allele in males has the same effect as a double allele dose in females. Default is "double", which corresponds to X-inactivation.} 
	\item{BR.girls}{To be used with \code{xchrom = TRUE}. Gives the ratio of baseline risk for females to the baseline risk for males.}
	\item{dire}{Gives the directory of the simulated data files. Missing by default, which means that none of the files are saved to files.}
	\item{ask}{Logical. If \code{TRUE}, hapSim will ask before overwriting the files in an already existing directory.}
	\item{cpus}{Allows parallel processing of its analyses. The \code{cpus} argument should preferably be set to the number of available CPUs. If set lower, it will save some capacity for other processes to run. Setting it too high should not cause any serious problems.}
	\item{slaveOutfile}{Character. If \code{slaveOutfile = ""} (default), output from all running cores will be printed in the standard R session window. Alternatively, the output can be saved to a file by specifying the file path and name.}
	\item{...}{Arguments to be used by \code{\link[Haplin]{haplin}}, \code{\link[Haplin]{haplinSlide}} or \code{\link[Haplin]{haplinStrat}}.}
}


\details{\code{hapRun} applies \code{\link[Haplin]{haplin}}, \code{\link[Haplin]{haplinSlide}} or \code{\link[Haplin]{haplinStrat}}
on each data file simulated by \code{\link[Haplin]{hapSim}}. 
It provides simulations on various family designs, i.e., triads, case-control, the hybrid design, and all intermediate designs.
The simulated files may accomodate fetal effects, maternal effects and/or parent-of-origin effects.
\code{hapRun} allows simulation of both autosomal and X-chromosome markers, assuming Hardy-Weinberg equilibrium.
It also enables simulation and calculation of gene-environment interaction effects.

Details on how to implement the arguments listed above are provided by \code{\link[Haplin]{hapSim}} and the Examples section below.
The stratum specific arguments may be simplified if the number of strata is one, or if the arguments are equal across all strata.

\code{\link[Haplin]{haplin}}, \code{\link[Haplin]{haplinStrat}} and \code{\link[Haplin]{haplinSlide}} will run with default values unless otherwise specified by \code{hapRun}.
For example, if \code{hapfunc = "haplin"}, \code{\link[Haplin]{haplin}} will use \code{response = "free"} unless \code{response = "mult"} is explicitly given as an argument in \code{hapRun}. 
Moreover, triads with missing data are only included in the haplin analysis if the argument \code{use.missing} equals \code{TRUE} (default in \code{hapRun}). Please confer \url{https://haplin.bitbucket.io/docu/Haplin_power.pdf} for further details and examples.
 
For information on the arguments to be passed on to \code{\link[Haplin]{haplin}}, \code{\link[Haplin]{haplinStrat}} and \code{\link[Haplin]{haplinSlide}}, please consult their help pages.

Note that \code{RR.mat} and \code{RRstar.mat} and \code{RRcm} and \code{RRcf} are required for \code{\link[Haplin]{hapSim}} to simulate maternal and parent-of-origin effects, respectively.
To calculate these effects, however, arguments \code{maternal = TRUE} and/or \code{poo = TRUE} must be specified.

\code{gen.missing.cases} and \code{gen.missing.controls} are flexible arguments. By default, both equal NULL, which means that no missing data are generated at random.
If the arguments are single numbers, missing data are generated at random with this proportion for all cases and/or controls. 
If the arguments are vectors of length equal to the number of loci, missing data are generated with the corresponding proportion for each locus.
The arguments can also be matrices with the number of rows equal to the number of loci and three columns. 
Each row corresponds to a locus, and the columns correspond to mothers, fathers and children, respectively. 
}

\value{If \code{hapfunc = "haplin"}, \code{hapRun} returns a dataframe consisting of results from running \code{\link[Haplin]{haplin}} on each simulated file.
The first two columns are:

	\item{sim.no  }{The name of the directory from which the results are calculated, i.e., the simulation number}
	\item{row.no  }{The row number within each simulation}
	
\code{\link[Haplin]{haptable}} gives detailed information of the full dataframe. 

If \code{hapfunc = "haplinSlide"}, \code{hapRun} returns a list of which each element contains the results from a single run of \code{\link[Haplin]{haplinSlide}}.
Consult \code{\link[Haplin]{suest}} for a thorough description of the output. Note, however, that \code{hapfunc = "haplinSlide"} is currently only implemented for diallelic markers, and the reference category is always chosen to be the first haplotype (see \code{\link[Haplin]{hapSim}} for a description of the haplotype grid). 

If \code{hapfunc = "haplinStrat"}, \code{\link[Haplin]{haplinStrat}} is used to estimate gene-effects in each stratum of the exposure covariate, and the results from all strata are compared using \code{\link[Haplin]{gxe}}.  \code{hapRun} returns a list, where each element is the result of a single run of 
\code{\link[Haplin]{gxe}}.

Additionaly, if \code{dire} is not missing by default, the simulated files from which the Haplin results are calculated, are stored in the given directory.}


\references{
Web Site: \url{https://haplin.bitbucket.io}
}

\author{Miriam Gjerdevik,\cr
with Hakon K. Gjessing\cr
Professor of Biostatistics\cr
Division of Epidemiology\cr
Norwegian Institute of Public Health\cr

\email{hakon.gjessing@uib.no}}

\seealso{\code{\link[Haplin]{haplin}}, \code{\link[Haplin]{haplinSlide}}, \code{\link[Haplin]{hapSim}}, \code{\link[Haplin]{haptable}}, \code{\link[Haplin]{suest}}, \code{\link[Haplin]{hapPower}}, \code{\link[Haplin]{hapPowerAsymp}}}

\examples{
\dontrun{
## Simulate Haplin results from 100 files using the multiplicative model in haplin. 
## The files consist of fetal effects at two diallelic markers,
## corresponding to haplo.freq = rep(0.25, 4), RR = c(2,1,1,1) 
## and RRstar = c(1,1,1,1). That is, the first allele has a doubled risk 
## relative to the rest. The data consists of a combination of 
## 100 case triads and 100 control triads with no missing data.
## No environmental factors are considered, i.e. the number of strata is one.
hapRun(nall = c(2,2), n.strata = 1, cases = c(mfc=100), controls = c(mfc=100),
haplo.freq = rep(0.25,4), RR = c(2,1,1,1), RRstar = c(1,1,1,1), 
hapfunc = "haplin", response = "mult", n.sim = 100, dire = "simfiles", ask = FALSE)

## Simulate power from 100 files applying haplinStrat. 
## The files consist of fetal and maternal effects at two diallelic markers.
## The data is simulated for 500 case triads and 200 control families in the first stratum,
## and 500 case triads and 500 control trids in the second.
## The fetal effects vary across strata,
## whereas the maternal effects are the same.
## One percent of the case triads are missing at random in the second stratum.
hapRun(nall = c(2,2), n.strata = 2, cases = c(mfc=500),
controls = list(c(mfc=200),c(mfc=500)), haplo.freq = rep(0.25,4), maternal = TRUE, 
RR = list(c(1.5,1,1,1),c(1,1,1,1)), RRstar = c(1,1,1,1),
RR.mat = c(1.5,1,1,1), RRstar.mat = c(1,1,1,1), 
gen.missing.cases = list(NULL,0.01), use.missing = TRUE, hapfunc = "haplinStrat", 
n.sim = 100, ask = FALSE)

## Simulate Haplin results from 100 files using haplin. 
## The files consist of fetal effects at one diallelic locus, 
## corresponding to haplo.freq = rep(0.5,2), RR = c(1.5,1) and RRstar = c(1,1).
## We have a combination of 100 case triads and 
## 100 control triads with no missing data. 
## No environmental effects are considered.
hapRun(nall = c(2), n.strata = 1, cases = c(mfc=100), controls = c(mfc=100),
haplo.freq = rep(0.5,2), RR = c(1.5,1), RRstar = c(1,1),
hapfunc = "haplin", n.sim = 100, dire = "simfiles", ask = FALSE)
}
}