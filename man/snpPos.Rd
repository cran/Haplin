\name{snpPos}
\alias{snpPos}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Finding the column numbers of SNP identifiers/SNP numbers in a ped file}
\description{
 Gives the column numbers of SNP identifiers or SNP numbers in a standard ped file, calculated from the SNP's positions in the corresponding map file. The column numbers are sorted in the same order as \code{snp.select}.
 These positions may be useful when extracting a selection of SNPs from a ped file.
}
\usage{
snpPos(snp.select, map.file, blank.lines.skip = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{snp.select}{A character vector of the SNP identifiers (RS codes) or a numeric vector of the SNP numbers.}
	\item{map.file}{A character string giving the name and path of the standard map file to be used.  See Details for a description of the standard map format.}
	\item{blank.lines.skip}{Logical. If "TRUE" (default), \code{snpPos} ignores blank lines in \code{map.file}.}
}
\details{
  To extract certain SNPs from a standard ped file, one has to know their positions in the ped file.
This can be obtained from the corresponding map file.\cr

The map file should look something like this:
\preformatted{
Chromosome SNP-identifier Base-pair-position
1               RS9629043             554636
1              RS12565286             711153 
1              RS12138618             740098
}
Alternatively, the map file could contain four columns. The column values should then be:
Chromosome, SNP-identifier, Genetic-distance, Base-pair-position.\cr
A header must be added to the map file if this does not already exist.\cr
\cr
The format of the corresponding ped file should be something like this:
\preformatted{
1104  1104-1  1104-2  1104-3  1  2  4  1  3  2
1104  1104-2       0       0  1  1  4  1  2  2
1104  1104-3       0       0  2  1  0  0  0  0
1105  1105-1  1105-2  1105-3  2  2  1  1  2  2
1105  1105-2       0       0  1  1  1  1  2  2
1105  1105-3       0       0  2  1  1  1  3  2
}
The column values are: Family id, Individual id, Father's id, Mother's id, Sex (1 = male, 2 = female, alternatively: 1 = male, 0 = female), and Case-control status (1 = controls, 2 = cases, alternatively: 0 = controls, 1 = cases).\cr
Column 7 and onwards contain the genotype data, with alleles in separate columns. A ``0'' is used to denote missing data.\cr

}
\note{The function does not check if the map file is formatted correctly or if the map and ped file have the same number of SNPs. The corresponding positions of the SNPs in the ped file may not be correct if the ped file has a different format from the given example.}

\value{
	A vector of the column numbers of the SNP identifiers/SNP numbers in the ped file, sorted in the same order as given in \code{snp.select}.
} 
\references{
Web Site: \url{http://www.uib.no/smis/gjessing/genetics/software/haplin/}}
\author{Miriam Gjerdevik,\cr
with Hakon K. Gjessing\cr
Professor of Biostatistics\cr
Division of Epidemiology\cr
Norwegian Institute of Public Health\cr
\email{hakon.gjessing@fhi.no}}

\seealso{\code{\link[Haplin]{convertPed}}, \code{\link[Haplin]{lineByLine}}}
\examples{

\dontrun{

# Find the column numbers of the SNP identifiers "RS9629043" and "RS12565286" in a standard ped file
snpPos(snp.select = c("RS9629043", "RS12565286"), map.file = "mygwas.map")
}

}


