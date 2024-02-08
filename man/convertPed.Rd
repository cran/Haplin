\name{convertPed}

\alias{convertPed}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Convert large ped files by creating unique IDs, converting allele coding and extracting a selection of SNPs}

\description{NOTE: This function is probably less useful now that GenABEL is no longer used by Haplin. The function is used to prepare a ped file for loading into GenABEL. However, GenABEL requires unique individual IDs in the file, not only unique within family. Furthermore, numeric allele coding 1,2,3,4 is not accepted. To fix this, \code{convertPed} can be run prior to running \code{prepPed}. This will create unique IDs and do the necessary allele recoding, and possibly also select and reorder SNPs. \code{convertPed} will also update the corresponding map file.}

\keyword{internal}

\usage{
convertPed(ped.infile, map.infile, ped.outfile, map.outfile, create.unique.id = FALSE,
convert, snp.select = NULL, choose.lines = NULL, col.sep = " ",
ask = TRUE, blank.lines.skip = TRUE, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
	\item{ped.infile}{A character string giving the name of the standard ped file to be modified. The name of the file is relative to the current working directory, unless the file name contains an absolute path.\cr  See Details for a description of the standard ped format.}
	\item{map.infile}{A character string giving the name and path of the to-be-modified standard map file. Optional if snp.select = NULL. A description of the standard map format is given in the Details section.}
	\item{ped.outfile}{A character string of the name and path of the converted ped file.}
	\item{map.outfile}{A character string giving the name and path of the modified map file.}
	\item{create.unique.id}{Logical. If "TRUE", the function creates a unique individual ID.}
	\item{convert}{No default. The option "ACGT_to_1234" recodes the SNP alleles from A,C,G,T to 1,2,3,4, whereas "1234_to_ACGT" converts from 1,2,3,4 to A,C,G,T. If "no_recode", no conversion occurs.}
	\item{snp.select}{A character vector of the SNP identifiers (RS codes) or a numeric vector of the SNP numbers to be extracted. Default is "NULL", which means that all SNPs are selected without reordering among the SNPs. The RS codes or SNP numbers may be listed in any order. Reordering among the selected SNPs will occur in the modified files corresponding to this listing.}
	\item{choose.lines}{A numeric vector of lines to be selected from the ped file. If "NULL" (default), all lines are selected.}
	\item{col.sep}{Specifies the separator that splits the columns in \code{ped.infile}. By default, col.sep = " " (space). To split at all types of space or blank characters, set col.sep = "[[:space:]]" or col.sep = "[[:blank:]]".}
	\item{ask}{Logical. Default is "TRUE". If set to "FALSE", an already existing outfile will be overwritten without asking.}
	\item{blank.lines.skip}{Logical. If "TRUE" (default), \code{convertPed} ignores blank lines in \code{ped.infile} and \code{map.infile}.}
	\item{verbose}{Logical. Default is "TRUE", which means that the line number is displayed for each iteration, i.e. each line read and modified, in addition to the first ten columns of the converted line.}
}

\details{
\code{convertPed} assumes a standard ped file as input. 
The format of the ped file should look something like this:
\preformatted{
1104  1  2  3  1  2  4  1  3  2  1  1
1104  2  0  0  1  1  4  1  2  2  4  1
1104  3  0  0  2  1  0  0  0  0  0  0
1105  1  2  3  2  2  1  1  2  2  4  1
1105  2  0  0  1  1  1  1  2  2  1  1
1105  3  0  0  2  1  1  1  3  2  4  4
}
The column values are: Family ID, Individual ID, Father's ID, Mother's ID, Sex (1 = male, 2 = female, alternatively: 1 = male, 0 = female), and Case-control status (1 = controls, 2 = cases, alternatively: 0 = controls, 1 = cases).\cr
Column 7 and onwards contain the genotype data, with alleles in separate columns, two columns representing one SNP. A ``0'' is used to denote missing data.\cr
\cr

The corresponding map file should look something like this:
\preformatted{
Chromosome SNP-identifier Base-pair-position
1               RS9629043             554636
1              RS12565286             711153 
1              RS12138618             740098
}
Alternatively, the map file could contain four columns. The column values should then be:
Chromosome, SNP-identifier, Genetic-distance, Base-pair-position.\cr
A header must be added to the map file if this does not already have one.\cr\cr

After creating unique individual IDs and recoding the SNP alleles from 1,2,3,4 to A,C,G,T (using \code{convertPed} with options \code{create.unique.id = TRUE} and \code{convert = "1234_to_ACGT"}), 
the ped file above should look like this:
\preformatted{
1104  1104_1  1104_2  1104_3  1  2  T  A  G  C  A  A
1104  1104_2       0       0  1  1  T  A  C  C  T  A
1104  1104_3       0       0  2  1  0  0  0  0  0  0
1105  1105_1  1105_2  1105_3  2  2  A  A  C  C  T  A
1105  1105_2       0       0  1  1  A  A  C  C  A  A
1105  1105_3       0       0  2  1  A  A  G  C  T  T
}
}

\note{The function does not check if the ped or map file is formatted correctly. For instance, if the alleles follows the generic A/B Illumina coding, \code{convertPed} may still be used to create unique individual IDs and extract a selection of SNPs. Using \code{convert = "ACGT_to_1234"} would however, result in nonsense.}

\value{
  There is no useful output; the objective of \code{convertPed} is the converted ped file and the modified map file. 
}

\references{
%Purcell S. PLINK v1.07. Web Site: \url{http://pngu.mgh.harvard.edu/purcell/plink/}\cr\cr
Web Site: \url{https://haplin.bitbucket.io}
}

\author{Miriam Gjerdevik,\cr
with Hakon K. Gjessing\cr
Professor of Biostatistics\cr
Division of Epidemiology\cr
Norwegian Institute of Public Health\cr
\email{hakon.gjessing@uib.no}}

\seealso{\code{\link[Haplin]{lineByLine}}, \code{Haplin:::lineConvert}, \code{\link[Haplin]{snpPos}}}

\examples{
\dontrun{

# Create unique individual IDs and recode SNP alleles from 1,2,3,4 to A,C,G,T
convertPed(ped.infile = "mygwas.ped", map.infile = "mygwas.map",
ped.outfile = "mygwas_modified.ped", map.outfile = "mygwas_modified.map",
create.unique.id = TRUE, convert = "1234_to_ACGT", ask = TRUE)

}
}



