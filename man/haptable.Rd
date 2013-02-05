\name{haptable}
\alias{haptable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create haplin table }
\description{
Create a comprehensive table of haplin output
}
\usage{
 haptable(object) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A \code{haplin} object, i.e. the result of running \code{haplin}.}
}
\details{
	\code{haptable} extracts the most important information from a \code{haplin} object to produce a summary table. The table can then be saved with, for instance, \code{write.table}, making the results easily accessible to other applications. You can also use \code{output} to produce and save the same results.
}
\value{
	A dataframe is returned, with the following columns:
	\item{Original }{Number of triads}
	\item{After.rem.NA }{	Number of triads after removal of missing}
	\item{After.rem.Mend..inc. }{	Number of triads after removal of Mendelian inconsistencies}
	\item{After.rem.rare.haplos }{	Number of triads after removal of rare haplotypes}
	\item{alleles }{ Alleles at each SNP}
	\item{counts }{ Frequency counts for alleles}
	\item{HWE.pv }{ P-value for deviation from HWE}
	\item{haplo }{ Haplotypes found during estimation}
	\item{pv.overall }{ Overall p-value}
	\item{hfreq.est. }{ Estimated haplotype frequencies}
	\item{hfreq.lower }{ Lower 95\% CI for estimated haplotype frequencies}
	\item{hfreq.upper }{ Upper 95\% CI for estimated haplotype frequencies}
	\item{reference }{ Reference method. If ref.cat is used, the reference category is labeled "ref"}
	\item{RR.est. }{ Estimated single dose relative risk}
	\item{RR.lower }{ Lower 95\% CI for estimated single dose relative risk}
	\item{RR.upper }{ Upper 95\% CI for estimated single dose relative risk}
	\item{RR.p.value }{ P-values for individual single dose estimates}
	\item{RRdd.est. }{ Estimated double dose relative risk}
	\item{RRdd.lower }{ Lower 95\% CI for estimated double dose relative risk}
	\item{RRdd.upper }{ Upper 95\% CI for estimated double dose relative risk}
	\item{RRdd.p.value }{ P-values for individual double dose estimates}
	\item{RRm.est. }{ Estimated single dose relative risk for maternal haplotype}
	\item{RRm.lower }{ Lower 95\% CI for estimated single dose relative risk for maternal haplotype}
	\item{RRm.upper }{ Upper 95\% CI for estimated single dose relative risk for maternal haplotype}
	\item{RRm.p.value }{ P-values for individual single dose estimates for maternal haplotype}
	\item{RRmdd.est. }{ Estimated double dose relative risk for maternal haplotype}
	\item{RRmdd.lower }{ Lower 95\% CI for estimated double dose relative risk for maternal haplotype}
	\item{RRmdd.upper }{ Upper 95\% CI for estimated double dose relative risk for maternal haplotype}
	\item{RRmdd.p.value }{ P-values for individual double dose estimates for maternal haplotype}
	Note that the results from the maternal haplotype risks will not be present unless \code{maternal = T} was used in the \code{haplin} call. 
}
\references{Web Site: \url{http://www.uib.no/smis/gjessing/genetics/software/haplin/}}
\author{Hakon K. Gjessing\cr
Professor of Biostatistics\cr
Division of Epidemiology\cr
Norwegian Institute of Public Health\cr
\email{hakon.gjessing@fhi.no}}
\note{Further information is found on the web page
}
\seealso{\code{\link{haplin}}, \code{\link{output}}}

\examples{
\dontrun{

# Produce a table containing the most important output from haplin:
res <- haplin("data.dat", use.missing = T, maternal = T)
haptable(res)



}}




