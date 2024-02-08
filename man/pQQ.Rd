\name{pQQ}
\alias{pQQ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{QQ-plot with confidence intervals for a vector of p-values}
\description{
 Produces a QQ-plot of p-values. The x-axis is -log10 of the expected p-values (under a null hypothesis of no effects), the y-axis is -log10 values of the actual p-values. A (pointwise) confidence interval can be drawn, and names of snps/genes corresponding to the most significant ones can be added.
}
\usage{
pQQ(pvals, nlabs = 6, conf = 0.95, lim, mark = 0.05, ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{pvals}{A vector of p-values. }
	\item{nlabs}{The number of (most significant) p-values to be labeled using \code{names(pvals)}.}
	\item{conf}{The confidence level of the pointwise confidence band. The default is 0.95. The confidence intervals are computed under the assumption of the p-values being drawn independently from a uniform [0,1] distribution. To leave out the confidence interval, set this to FALSE.}
	\item{lim}{A vector of length 2 giving the plot limits (on a log10 scal, for instance c(0,4)). Plot limits are computed automatically. However, if other plot limits are desireable, they can be set using this argument. }
	\item{mark}{By default, the 0.05 significance level is marked by lines. Can be changed to a different value, or set to FALSE. }
	\item{...}{Other arguments passed on to the plotting function.}
}
\details{
  The \code{pvals} argument should be a vector of p-values to be plotted. If the vector has names corresponding to marker (snp) names, the plot will label some of the most significant points with the marker names.
}
\value{
	No value is returned.
}
\references{ Gjessing HK and Lie RT. Case-parent triads: Estimating single- and double-dose effects of fetal and maternal disease gene haplotypes. Annals of Human Genetics (2006) 70, pp. 382-396.\cr\cr
Web Site: \url{https://haplin.bitbucket.io}}
\author{Hakon K. Gjessing\cr
Professor of Biostatistics\cr
Division of Epidemiology\cr
Norwegian Institute of Public Health\cr
\email{hakon.gjessing@uib.no}}

