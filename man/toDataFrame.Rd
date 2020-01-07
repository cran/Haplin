\name{toDataFrame}
\alias{toDataFrame}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Stack dataframes from haplinSlide into a single dataframe}
\description{
When haplinSlide is run with the option \code{table.output = T}, the result is a list of \code{haptable}s, i.e. tables with summary haplin results for each window \code{haplinSlide} is run on. \code{toDataFrame} stacks the separate dataframes into one large dataframe containing all results.
}
\usage{
toDataFrame(x, reduce = F) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{x}{The output from \code{haplinSlide} run with option \code{table.output = TRUE}.}
	\item{reduce}{Reduce output to one line per marker}
}
\details{
When \code{haplinSlide} is run with \code{winlength = 1} on SNP markers, each table in the output has only two rows, and can be condensed to a single row. By setting the argument \code{reduce} to TRUE, \code{toDataFrame} reduces each table to one line and returns a dataframe with one line for each SNP. In more general situations, with multi-allelic markers or, more commonly, \code{winlength} set to 2 or more, each output table will typically have more than two rows and cannot be reduced, so \code{reduce} should be set to FALSE.
}
\value{
	The output is a dataframe. First column contains the marker names. Second column are row numbers, counted within each output table. The remaining columns are identical to the individual output columns, which are described in more detail in the help file for \code{\link{haptable}}.
}
\references{ Gjessing HK and Lie RT. Case-parent triads: Estimating single- and double-dose effects of fetal and maternal disease gene haplotypes. Annals of Human Genetics (2006) 70, pp. 382-396.\cr\cr
Web Site: \url{http://folk.uib.no/gjessing/genetics/software/haplin/}}
\author{Hakon K. Gjessing\cr
Professor of Biostatistics\cr
Division of Epidemiology\cr
Norwegian Institute of Public Health\cr
\email{hakon.gjessing@uib.no}}

