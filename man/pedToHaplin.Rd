\name{pedToHaplin}
\alias{pedToHaplin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert from ped format data to Haplin format}
\description{
 Converts an ASCII file from a standard ped format to the Haplin format
}
\usage{
pedToHaplin(indata, outdata, merge = F, na.strings = "0", sep, 
colnames.out = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{indata}{A character string giving the name and path of the ASCII data file to be converted. }
	\item{outdata}{A character string giving the name and path for saving the converted file. }
	\item{merge}{If the alleles of each genotype are in two separate columns in the \code{indata} file, they must be merged (with the ";" separator) in the \code{outdata} file. This is done by setting \code{merge = TRUE}. Otherwise, it must be set to \code{FALSE}.}
	\item{na.strings}{The symbol used to denote missing data in \code{indata}. It is passed directly to R's \code{read.table}}
	\item{sep}{Column separator in \code{indata}. If unspecified, any white space will be used, as in \code{read.table}.} 
	\item{colnames.out}{Provided just for the purpose of checking data. If TRUE, adds colnames to the \code{outdata} file to make it more readable. NOTE: Haplin does currently not use colnames, so this should be set to FALSE when producing the file to run on.}
}
\details{
  Important: The first 6 columns should always be family id, individual id, father's id, mother's id, sex and casetype, in that order, then followed by the genetic data columns. If the genetic data columns are separated into two individual alleles, one should use the option \code{merge = TRUE} to merge them in the output file. If they are already joined in single columns, for instance as CT or C;T, \code{merge} should be set to \code{FALSE} (default).

Additional covariates can be included in the input file. If so, they should be placed after the 6 leading columns but before the genetic data. In this case, one should make sure the genetic data columns are already merged, and that \code{merge = FALSE}. (The \code{merge = TRUE} option when covariates are present will hopefully be implemented at some point...)

Note that the output file usually has three columns before (to the left of) the columns containing genetic data. These columns are family id, sex, and casetype. When running haplin on the output file one should specify the argument 'n.vars = 3' in haplin. If the data are from the x chromosome the haplin arguments should also include 'sex = 2' and 'xchrom = T'. Similarly, if the casetype variable is a case-control indicator one should use the argument 'ccvar = 3'. If the intention is to only run haplin on the cases the case triads should be saved separately in a new file prior to running haplin on it.
}
\value{
	The \code{outdata} file is written to disk. \code{pedToHaplin} returns (invisibly) the converted data file.
}
\references{ Gjessing HK and Lie RT. Case-parent triads: Estimating single- and double-dose effects of fetal and maternal disease gene haplotypes. Annals of Human Genetics (2006) 70, pp. 382-396.\cr\cr
Web Site: \url{https://people.uib.no/gjessing/genetics/software/haplin/}}
\author{Hakon K. Gjessing\cr
Professor of Biostatistics\cr
Division of Epidemiology\cr
Norwegian Institute of Public Health\cr
\email{hakon.gjessing@uib.no}}

\section{Warning}{Data files come in many shapes and formats, so you should always check the output from \code{pedToHaplin} before using it.}

\note{
	Further information is found on the web page.
}
\seealso{\code{\link{haplin}}}
\examples{

\dontrun{

# Standard run on supplied test file:
pedToHaplin("test_pedToHaplin.ped", outdata = "test_pedToHaplin_result.txt", 
colnames.out = F, merge = T)

}

}
