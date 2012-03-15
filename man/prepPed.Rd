\name{prepPed}
\alias{prepPed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract family and phenotype information from a ped-format file, to prepare for use in Haplin}
\description{
 Creates a pedIndex file containing family information, a phenotype file, and optionally an ``dummy'' map file. The files are used by GenABEL when loading data into R, and by Haplin when converting from a GenABEL file to a Haplin file. 
}
\usage{
prepPed(pedfile, outdir, create.map = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{pedfile}{A character string giving the name and path of the ped-format file to be used. }
	\item{outdir}{The directory where the pedIndex file, phenotype file, and optionally the map file should be saved.}
	\item{create.map}{Logical. If "TRUE", \code{prepPed} creates a dummy map file which can be used by GenABEL when loading data into R. Can be used if no map file is available.}
}
\details{
  To use Haplin on a large ped-format file, it should first be converted to a GenABEL raw file and loaded into R. Since GenABEL does not retain family information available in the ped file, \code{prepPed} should first be run on the file to extract the necessary family and phenotype information. \code{prepPed} stores family information in a .pedIndex file with the same name as the ped file, and saves it in the \code{outdir} directory. Similarly, it creates a phenotype file (.ph), which contains the individual ID, the sex variable, and the case-control status. Optionally, it can construct a simple .map file, which can be used in situations where no real map file (corresponding to the ped file) is available.\cr

The format of the ped file should be something like this:
\preformatted{
1104  1104-1  1104-2  1104-3  1  1  A  B  B  B
1104  1104-2       0       0  1  1  B  B  A  B
1104  1104-3       0       0  2  1  A  B  A  B
1105  1105-1  1105-2  1105-3  2  2  B  B  A  A
1105  1105-2       0       0  1  2  B  B  A  A
1105  1105-3       0       0  2  2  0  0  A  A
}
The column values are: Family id, Individual id, Father's id, Mother's id, Sex (1 = male, 2 = female), and Case-control status (0 = controls, 1 = cases).\cr

Column 7 and onwards contain the genotype data, with alleles in separate columns, or joined, as AB BB, etc. A ``0'' is used to denote missing data.\cr


Missing values in the sex and case-control columns are not accepted.
}
\note{More details on input format, output format etc. is found in the Haplin data description on the web page.}

\value{
  There is no useful output; the task of \code{prepPed} is to save the extracted information in the \code{outdir} directory.
}
\references{ Gjessing HK and Lie RT. Case-parent triads: Estimating single- and double-dose effects of fetal and maternal disease gene haplotypes. Annals of Human Genetics (2006) 70, pp. 382-396.\cr\cr
Web Site: \url{http://www.uib.no/smis/gjessing/genetics/software/haplin/}}
\author{Hakon K. Gjessing\cr
Professor of Biostatistics\cr
Division of Epidemiology\cr
Norwegian Institute of Public Health\cr
\email{hakon.gjessing@fhi.no}}

\seealso{\code{\link[GenABEL]{convert.snp.ped}}, \code{\link[GenABEL]{load.gwaa.data}}}
\examples{

\dontrun{

# Create the files mygwas.pedIndex, mygwas.ph and mygwas.map in the "data" directory
prepPed(pedfile = "data/mygwas.ped", outdir = "data", create.map = T)

}

}
