\name{gwaaToHaplin}
\alias{gwaaToHaplin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert a loaded GenABEL \code{gwaa.data} object to a Haplin format character matrix}
\description{
 Mostly for internal use. Can be useful as a check that a loaded GenABEL data object converts corretly to a Haplin format; otherwise not needed by the user. 
}
\usage{
  gwaaToHaplin(data, pedIndex, design = "triad")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{data}{Same as \code{data} argument in \code{haplin}.}
	\item{pedIndex}{Same as \code{pedIndex} argument in \code{haplin}.}
	\item{design}{Same as \code{design} argument in \code{haplin}.}
}
\details{
  To use Haplin on a large ped-format file, it should first be converted to a GenABEL raw file and loaded into R. See the documentation for \code{haplin} how this is done. Once a GenABEL object has been loaded into R, it can be fed directly to Haplin, and the internal conversion to Haplin format is done by \code{gwaaToHaplin}. The user might want to check the conversion and can apply \code{gwaaToHaplin} to the loaded data object to see the Haplin data matrix. Apart from that, the user has little need for \code{gwaaToHaplin}.
}
\section{Warning}{Do not try to convert a large GenABEL object directly. You should rather select a small subset of markers for inspection.}

\value{
  A data matrix (character) with one row for each family. The column names should be more or less self-explanatory. "m", "f", and "c" are used to names columns belonging to the mother, the father, and the child, respectively. Alleles are in separate columns.
}
\references{ Gjessing HK and Lie RT. Case-parent triads: Estimating single- and double-dose effects of fetal and maternal disease gene haplotypes. Annals of Human Genetics (2006) 70, pp. 382-396.\cr\cr
Web Site: \url{http://folk.uib.no/gjessing/genetics/software/haplin/}}
\author{Hakon K. Gjessing\cr
Professor of Biostatistics\cr
Division of Epidemiology\cr
Norwegian Institute of Public Health\cr
\email{hakon.gjessing@fhi.no}}

\seealso{\code{\link{haplin}}, \code{\link{prepPed}}, \code{\link[GenABEL]{convert.snp.ped}}, \code{\link[GenABEL]{load.gwaa.data}}}
\examples{

\dontrun{

## Extract family and phenotype information:
prepPed(pedfile = "data/mygwas.ped", outdir = "data", create.map = T)

## Convert to raw file format:
convert.snp.ped(pedfile = "data/mygwas.ped", 
 mapfile = "data/mygwas.map", outfile = "data/mygwas.raw")

## Load into R:
mygwas.data <- load.gwaa.data(phenofile = "data/mygwas.ph", 
 genofile = "data/mygwas.raw")

## Convert subset of data:
temp <- gwaaToHaplin(mygwas.data[, 1:3], 
 pedIndex = "data/mygwas.pedIndex")

## Inspect:
head(temp)


}

}
