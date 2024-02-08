\name{cbindFiles}
\alias{cbindFiles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Combine a sequence of files by columns}
\description{
Takes a sequence of files and combines them column-wise (side-by-side), i.e. reads each file line by line, pastes corresponding lines, then writes to outfile. Combining files line by line is especially useful when working with large datasets, where the reading of entire files may be time consuming and require a large amount of memory.
}
\usage{
cbindFiles(infiles, outfile, col.sep, ask = TRUE, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{infiles}{A character vector of names (and paths) of the files to combine.}
	\item{outfile}{A character string giving the name of the resulting file. The name of the file is relative to the current working directory, unless the file name contains a definite path.}
	\item{col.sep}{Specifies the column separator which will be inserted between files.}
	\item{ask}{Logical. Default is "TRUE". If set to "FALSE", an already existing outfile will be overwritten without asking.}
	\item{verbose}{Logical. Default is "TRUE", which means that the line number is displayed for each iteration, i.e. each combined line.}
}
\details{
	The function \code{\link[base]{cbind}} combines R objects by columns. However, reading large datafiles may require a large amount of memory and be extremely time consuming.
	\code{cbindFiles} combines the files column-wise, one line at a time, and then writes each line to \code{outfile}. This avoids reading the full files into memory.\cr
	If \code{infiles} contains only one file, this file will be copied directly to \code{outfile}, i.e. without any line-by-line modifications.
}
\value{
There is no useful output; the objective of \code{cbindFiles} is to produce \code{outfile}.
}
\references{Web Site: \url{https://haplin.bitbucket.io}}
\author{
Miriam Gjerdevik,\cr
with Hakon K. Gjessing\cr
Professor of Biostatistics\cr
Division of Epidemiology\cr
Norwegian Institute of Public Health\cr
\email{hakon.gjessing@uib.no}
}

\seealso{\code{\link[Haplin]{rbindFiles}}, \code{\link[Haplin]{lineByLine}}}

\examples{
\dontrun{

# Combines the three infiles side-by-side
cbindFiles(infiles = c("myfile1.txt", "myfile2.txt", 
"myfile3.txt"), outfile = "myfile_combined_by_columns.txt", 
col.sep = " ", ask = TRUE, verbose = TRUE)

}
}

