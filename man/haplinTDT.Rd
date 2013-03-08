\name{haplinTDT}
\alias{haplinTDT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Transmission disequilibrium tests for case-parent triad data}
\description{
  Produces an object of class \code{haplinTDT}, which is the result of
  Transmission disequilibrium tests of data.
}
\usage{
haplinTDT(filename, nsim.perm = 0, select.gender = NULL, 
method = c("tdt", "hhrr", "trimm"), names.marker = NULL, 
use.haplotypes = FALSE, use.ambiguous = TRUE, design = "triad", 
markers = "ALL", n.vars = 0, sep = " ", allele.sep = ";", 
na.strings = "NA", use.missing = FALSE, xchrom = FALSE, 
sex = NULL, threshold = 0.01, verbose = TRUE, printout = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{}{Of the following arguments, only \code{filename} is required. Use of the remaining arguments will depend on the type of analysis. }
  \item{filename}{A character string giving the name and path of
    the ASCII data file to be read. }
  \item{nsim.perm}{Number of permutations. Default is 0, which means that \code{haplinTDT} does not do
    a permutation test.}
  \item{select.gender}{Do the analysis for a gender subset. Values: 1, 2, or NULL.
    1: Male, 2: Female, NULL: All. Default is NULL.}
  \item{method}{A character vector containing the methods that are used
    int the analysis. Possible values are "tdt", "hhrr" and
    "trimm". Default are all three.}
  \item{names.marker}{Marker names. Default is NULL which means that the
    markers are denoted 1, 2, ..., # markers.}
  \item{use.haplotypes}{A logical value, default is FALSE. If
    use.haplotypes=TRUE haplotypes corresponding to the individual
    markers are reconstructed by \code{\link{haplin}} using the EM
    algorithm.
    The haplotypes are then
    analysed as a single multiallelic marker. If use.haplotypes=FALSE,
    the markers are analysed individually.}
  \item{use.ambiguous}{A logical value, default is TRUE. If FALSE then
    we remove those family triads where it is ambiguous which allele is
    transferred to the child from a parent.}
  \item{design}{For the moment only the value "triad" is allowed. It is used for the standard case triad
    design, without independent controls.}
  \item{markers}{Default is "ALL", which means \code{haplinTDT} uses all
    available markers in the data set in the analysis. If
    \code{use.haplotypes = TRUE} then for the current
    version of \code{\link{haplin}} the number of markers used at a single run
    should probably not exceed 4 or 5 due to the computational
    burden. The markers argument can be used to select appropriate
    markers from the file without creating a new file for the selected
    markers. For instance, if markers is set to c(2,4), \code{haplinTDT} will
    only use the second and fourth markers supplied in the data
    set. When running \code{haplinTDT}, it may be a good idea to start
    exploring a few markers at a time, using this argument.  }
  \item{n.vars}{Numeric. The number of variables (columns) in the data
    file before (to the left) of the genetic data.}
  \item{sep}{The character separator used in the data file to separate
    between "columns", where each column contains the two alleles of a
    single individual at a single marker. }
  \item{allele.sep}{The character separator used in the data file to
    separate the two alleles for a single individual in a single
    marker. The recommended (default) separator is ";", but for SNPs an
    empty "" is also common.}
  \item{na.strings}{The character string indicating missing data in the
    data file. Default is to use "NA" in place of, for instance, C;T for
    a SNP that hasn't been typed in that individual.}
  \item{use.missing}{A logical value used to determine whether triads
    with missing data should be included in the analysis. When set to
    TRUE, \code{\link{haplinTDT}}  uses \code{\link{haplin}} to reconstruct the markers or haplotypes. The default, however, is
    FALSE. When FALSE, all triads having  any sort of missing data are
    excluded before the analysis is run. Note that  \code{haplinTDT} only looks at
    markers actually used in the analysis,  so that if the markers
    argument (see below) is used to select a collection of markers for
    analysis, \code{haplinTDT} only excludes triads with missing data on the
    included markers.}
  \item{xchrom}{Logical, defaults to "FALSE". If set to "TRUE",
    \code{haplinTDT} assumes the markers are on the x-chromosome. This option
    should be combined with specifying the \code{sex} argument.}
  \item{sex}{A numeric value specifying which of the data columns that
    contains the sex variable. The variable should be coded 1 for males
    and 2 for females. To be used with \code{xchrom = TRUE}.}
  \item{threshold}{Sets the (approximate) lower limit for the haplotype
    frequencies of those haplotypes that should be retained in the
    analysis. Hapotypes that are less frequent are removed, and
    information about this is given in the output. }
  \item{verbose}{Default is T (=TRUE). }%During the EM algorithm, haplinTDT prints the estimated parameters and deviance for each step. To avoid the output, set this argument to F (=FALSE). }
  \item{printout}{Logical. If TRUE (default), \code{haplinTDT} prints a full
    summary of the results after finishing the estimation. If FALSE, no
    such printout is given, but the \code{summary} function can later be
    applied to a saved result to get the same summary.}
%  \item{...}{Remaining arguments to be used by\code{\link{haplin}} in each run.}
}
\details{
  Threee types of transmission disequilibrium tests (TDT) are
  provided:

  Let t_ij be the number of parents transmitting
  allele i and not j to its child and let n be number of alleles. The
  standard TDT test is then defined as
  the sum of terms (t_ij - t_ji)^2/(t_ij + t_ji) for 1<=i<j<=n.
  This sum is asymptotically chi-squared with n(n-1) degrees of freedom when
  the marker and disease loci are unlinked or not associated.
  
  Let t_i. and t_.i  be the marginal
  totals of t_ij. The Haplotype-based  haplotype relative risk (HHRR) test is
  then defined as
  the sum of (t_i. - t_.i)^2/(t_i. + t_.i) for i = 1, 2, ..., n.
  The HHRR test statistic is asymptotically chi-squared with n-1 degrees of freedom.

  The Triad Multi-Marker test (TRIMM) test is only defined for
  diallelic markers.

  If use.ambiguous = FALSE, then all ambiguous trios will be removed.
  Otherwise, the different contributions
  to TDT, HHRR and TRIMM are weighted with the probabilities of the
  different transmission configurations of alleles from parent to child. For example if the
  parents and the child are all heterozygous 1/2, then with probability
  0.5 the mother (or father) will transfer allele 1 and not allele 2.
  The standard formulation of the TDT and HHRR tests
  correspond to having use.ambiguous = TRUE.
}

\value{
  An object of class \code{haplinTDT} is returned
}
\references{Haplotype-based haplotype haplotype relative risk (HHRR):\cr
  Terwilliger JD and Ott J. A haplotype-based 'haplotype relative risk'
  approach to detecting allelic associations. Human Heredity (1992)
  42(6), pp. 337-46.\cr

  Transmission disequilibrium test (TDT):\cr
  Spielman RS, McGinnis RE and Ewens WJ. Transmission test for linkage
  disequilibrium: the insulin gene region and insulin-dependent diabetes
  mellitus (IDDM). American Journal of Human Genetics (1993) 52(3),
  pp. 506-16.\cr

  Triad multi-marker test (TRIMM):\cr
  Shi M, Umbach DM and Weinberg CR. Identification of
  risk-related haplotypes with the use of multiple SNPs from nuclear
  families. The American Journal of Human Genetics (2007) 81,
  pp. 53-66.
  }

\author{Oivind Skare,\cr with Hakon K. Gjessing\cr
  Division of Epidemiology\cr
  Norwegian Institute of Public Health\cr
  \email{hakon.gjessing@fhi.no}}

%\section{Warning}{Some warnings????.}

\note{Further information is found on the web page.
}
\seealso{\code{\link{summary.haplinTDT}}, \code{\link{plot.haplinTDT}}}
\examples{

\dontrun{

# Standard run with permutation test:
res <- haplinTDT("data.dat", nsim.perm=1000)
# Plot the saved result:
plot(res)
# A full summary of saved result including p-values
summary(res)

# Include missing values:
res <- haplinTDT("data.dat", nsim.perm=1000, use.missing=TRUE)
}

}
