\name{rbindFiles}
\alias{rbindFiles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Combine a sequence of files by rows}

\description{Takes a sequence of files and combines them by rows, without reading the full files into memory. This is especially useful when dealing with large datasets, where the reading of entire files may be time consuming and require a large amount of memory.}

\usage{
rbindFiles(infiles, outfile, col.sep, header = FALSE, ask = TRUE, 
verbose = FALSE, add.file.number = FALSE, blank.lines.skip = FALSE)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
	\item{infiles}{A character vector of names (and paths) of the files to combine.}
	\item{outfile}{A character string giving the name of the modified file. The name of the file is relative to the current working directory, unless the file name contains a definite path.}
	\item{col.sep}{Specifies the separator used to split the columns in the files. To split at all types of spaces or blank characters, set \code{col.sep = "[[:space:]]"} or \code{col.sep = "[[:blank:]]"}.}
	\item{header}{A logical variable which indicates if the first line in each file contains the names of the variables. If "TRUE", \code{outfile} will display this header in its first row, assuming the headers for each file are identical. Equals FALSE by default, i.e. no headers assumed.}
	\item{ask}{Logical. Default is "TRUE". If set to "FALSE", an already existing outfile will be overwritten without asking.}
	\item{verbose}{Logical. Default is "TRUE", which means that the line number is displayed for each iteration, i.e. each combined line.}
	\item{add.file.number}{A logical variable which equals "FALSE" by deafult. If "TRUE", an extra first column will be added to the outfile, consisting of the file numbers for each line.}
	\item{blank.lines.skip}{Logical. If "TRUE" (default), \code{\link[Haplin]{lineByLine}} ignores blank lines in the input.}
}

\details{
	The function \code{\link[base]{rbind}} combines R objects by rows. However, reading large data files may require a large amount of memory and be extremely time consuming.
	\code{rbindFiles} avoids reading the full files into memory. It reads the files line by line, possibly modifies each line, then writes to outfile. 
	If however, \code{header}, \code{verbose}, \code{add.file.number} and \code{blank.lines.skip} are all set to "FALSE" 
	(their default values), the files are appended directly, thus evading line-by-line modifications. 
	In the case where \code{infiles} contains only one file and no output or modifications are requested 
	(\code{verbose}, \code{add.file.number} and \code{blank.lines.skip} equal "FALSE"), an identical copy of this file is made.
}

\note{Combining the files by reading each file line by line is less time efficient than appending the files directly. For this reason, if \code{header = FALSE}, changing the values of the logical variables \code{verbose}, \code{add.file.number}
and \code{blank.lines.skip} from "FALSE" to "TRUE" should not be done unless absolutely necessary.
}

\value{There is no useful output; the objective of \code{rbindFiles} is to produce \code{outfile}.}

\references{
Web Site: \url{http://www.uib.no/smis/gjessing/genetics/software/haplin/}
}

\author{Miriam Gjerdevik,\cr
with Hakon K. Gjessing\cr
Professor of Biostatistics\cr
Division of Epidemiology\cr
Norwegian Institute of Public Health\cr
\email{hakon.gjessing@fhi.no}}

\seealso{\code{\link[Haplin]{cbindFiles}}, \code{\link[Haplin]{lineByLine}}}

\examples{
\dontrun{

# Combines the three infiles, by rows
rbindFiles(file.names = c("myfile1.txt", "myfile2.txt", "myfile3.txt"), 
outfile = "myfile_combined_by_rows.txt", col.sep = " ", header = TRUE, verbose = TRUE)

}
}
