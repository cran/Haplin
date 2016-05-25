\name{hapPower}

\alias{hapPower}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Power simulation for association analyses with Haplin}

\description{Simulates the statistical power of genetic analyses assessing fetal effects, maternal effects and/or parent-of-origin effects. Effects of X-chromosome genes and gene-environment interaction effects are also allowed. 
}

\usage{
hapPower(hapRun.result, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
	\item{hapRun.result}{The result of running \code{\link[Haplin]{hapRun}}}
	\item{alpha}{\code{alpha} is the Type I Error probablity. Equals 0.05 by default.}
}


\details{The Haplin framework includes different modules for assessing genetic effects: \code{\link[Haplin]{haplin}}, 
\code{\link[Haplin]{haplinStrat}} and \code{\link[Haplin]{haplinSlide}}. \code{hapPower} simulates the power of these analyses, which enables power calculations of fetal effects, maternal effects and/or parent-of-origin effects. Various family designs, i.e., triads, case-control, the hybrid design, and all intermediate designs, are possible.
It also allows power calculation of gene-environment interaction effects and effects on X-chromosome markers.

\code{hapPower} calculates statistical power using the result of \code{\link[Haplin]{hapRun}}, and the target effects must be specified in this function, see Examples below, and details in \url{http://folk.uib.no/gjessing/genetics/software/haplin/Haplin_power.pdf}.
}

\value{\code{hapPower} returns the simulated power.}

\references{
Web Site: \url{http://folk.uib.no/gjessing/genetics/software/haplin/}
}

\author{Miriam Gjerdevik,\cr
with Hakon K. Gjessing\cr
Professor of Biostatistics\cr
Division of Epidemiology\cr
Norwegian Institute of Public Health\cr

\email{hakon.gjessing@fhi.no}}

\seealso{\code{\link[Haplin]{haplin}}, \code{\link[Haplin]{haplinSlide}}, \code{\link[Haplin]{hapSim}}, \code{\link[Haplin]{hapRun}}, \code{\link[Haplin]{snpPower}}}

\examples{

\dontrun{

## Simulate power from 100 files using haplin. 
## The files consist of fetal effects at two diallelic markers,
## corresponding to haplo.freq = rep(0.25, 4), RR = c(2,1,1,1) and RRstar = c(1,1,1,1).
## The power is simulated for the combination of 100 case triads 
## and 100 control triads with no missing data at a 0.05 significance level,
## applying a multiplicative model.
.hapRun <- hapRun(nall = c(2,2), n.strata = 1, cases = c(mfc=100), controls = c(mfc=100),
haplo.freq = rep(0.25,4), RR = c(2,1,1,1), RRstar = c(1,1,1,1), 
hapfunc = "haplin", response = "mult", n.sim = 100, dire = "simfiles")
hapPower(.hapRun)

## Simulate power from 100 files applying haplinStrat. 
## The files consist of fetal and maternal effects at two diallelic markers.
## The data is simulated for 500 case triads and 200 control families in the first stratum,
## and 500 case triads and 500 control trids in the second.
## The fetal effects vary across strata,
## whereas the maternal effects are the same.
## One percent of the case triads are missing at random in the second stratum.
.hapRun <- hapRun(nall = c(2,2), n.strata = 2, cases = c(mfc=500),
controls = list(c(mfc=200),c(mfc=500)), haplo.freq = rep(0.25,4), maternal = T, 
RR = list(c(1.5,1,1,1),c(1,1,1,1)), RRstar = c(1,1,1,1),
RR.mat = c(1.5,1,1,1), RRstar.mat = c(1,1,1,1), gen.missing.cases = list(NULL,0.01), 
use.missing = T, hapfunc = "haplinStrat", n.sim = 100)
hapPower(.hapRun)

## Simulate power at the 10% significance level from 100 files using haplin.
## The files consist of fetal effects at one diallelic locus,
## corresponding to haplo.freq = c(0.1,0.9), RR = c(1.5,1) and RRstar = c(1,1).
## The data consists of a combination of 100 case triads and 100 control triads.
.hapRun <- hapRun(nall = c(2), cases = c(mfc=100), controls = c(mfc=100),
haplo.freq = c(0.1,0.9), RR = c(2,1), RRstar = c(1,1),
hapfunc = "haplin", response = "mult", n.sim = 1000)
hapPower(.hapRun, alpha= 0.10)

## The latter example, applying response = "mult", should be comparable to 
## the theoretic calculations of snpPower.
snpPower(cases = list(mfc=100), controls = list(mfc=100),
RR = 2, MAF = 0.1, alpha = 0.10)
}
}
