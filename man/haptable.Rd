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
	----
	A dataframe is returned, with the following columns:
	----
	\item{marker}{Name(s) of marker(s) investigated}
	\item{alleles}{A listing of the alleles found at each marker}
	\item{counts}{Frequency counts of alleles at each marker}
	\item{HWE.pv}{P-value from Hardy-Weinberg equilibrium test at each marker}
	\item{Original}{Number of triads before removal}
	\item{After.rem.NA}{Number of triads after removal of missing}
	\item{After.rem.Mend.inc.}{Number of triads after removal of Mendelian inconsistencies}
	\item{After.rem.unused.haplos}{Number of triads after removal of unused (rare) haplotypes}
	\item{pv.overall}{Overall likelihood ratio p-value (test of all genetic effects combined)}
	\item{haplos}{Haplotypes (or single-marker alleles) found during estimation}
	\item{haplofreq}{Estimated haplotype frequencies}
	\item{haplofreq.lower}{Lower 95\% CI for estimated haplotype frequencies}
	\item{haplofreq.upper}{Upper 95\% CI for estimated haplotype frequencies}
	\item{reference}{Reference method. If ref.cat is used, the reference category is labeled "ref"}
	\item{RR.est.}{Estimated single dose relative risk}
	\item{RR.lower}{Lower 95\% CI for single dose relative risk}
	\item{RR.upper}{Upper 95\% CI for single dose relative risk}
	\item{RR.p.value}{P-values for individual single dose effect}
	\item{RRdd.est.}{Estimated double dose relative risk}
	\item{RRdd.lower}{Lower 95\% CI for double dose relative risk}
	\item{RRdd.upper}{Upper 95\% CI for double dose relative risk}
	\item{RRdd.p.value}{P-values for individual double dose effect}
	\item{NOTE1}{When maternal = TRUE, there will be additional columns:}
	\item{RRm.est.}{Estimated single dose relative risk for maternal haplotype}
	\item{RRm.lower }{Lower 95\% CI for single dose relative risk for maternal haplotype}
	\item{RRm.upper}{Upper 95\% CI for single dose relative risk for maternal haplotype}
	\item{RRm.p.value}{P-values for individual single dose effect of maternal haplotype}
	\item{RRmdd.est.}{Estimated double dose relative risk for maternal haplotype}
	\item{RRmdd.lower}{Lower 95\% CI for double dose relative risk for maternal haplotype}
	\item{RRmdd.upper}{Upper 95\% CI for double dose relative risk for maternal haplotype}
	\item{RRmdd.p.value}{P-values for individual double dose effect of maternal haplotype}
	\item{NOTE2}{When poo = TRUE, the RR.est., RR.lower, and RR.upper columns will be replaced by the following columns:}
	\item{RRcm.est.}{Estimated single dose relative risk, when inherited from the mother}
	\item{RRcm.lower}{Lower 95\% CI for single dose relative risk, when inherited from the mother}
	\item{RRcm.upper}{Upper 95\% CI for single dose relative risk, when inherited from the mother}
	\item{RRcm.p.value}{P-values for individual single dose effects, when inherited from the mother}
	\item{RRcf.est.}{Estimated single dose relative risk, when inherited from the father}
	\item{RRcf.lower}{Lower 95\% CI for single dose relative risk, when inherited from the father}
	\item{RRcf.upper}{Upper 95\% CI for single dose relative risk, when inherited from the father}
	\item{RRcf.p.value}{P-values for individual single dose effects, when inherited from the father}
	\item{RRcm_RRcf.est.}{An estimate of parent-of-origin effect, i.e. the ratio RRcm/RRcf}
	\item{RRcm_RRcf.lower}{Lower 95\% CI for ratio RRcm/RRcf}
	\item{RRcm_RRcf.upper}{Upper 95\% CI for ratio RRcm/RRcf}
	\item{RRcm_RRcf.p.value}{P-value for parent-of-origin effect RRcm/RRcf at that marker}	   
}
\references{Web Site: \url{https://haplin.bitbucket.io}}
\author{Hakon K. Gjessing\cr
Professor of Biostatistics\cr
Division of Epidemiology\cr
Norwegian Institute of Public Health\cr
\email{hakon.gjessing@uib.no}}
\note{Further information is found on the web page
}
\seealso{\code{\link{haplin}}, \code{\link{output}}}
\examples{
\dontrun{

# Produce a table containing the most important output from haplin:
res <- haplin("data.dat", use.missing = T, maternal = T)
haptable(res)
}}
