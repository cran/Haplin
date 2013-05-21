\name{lineByLine}
\alias{lineByLine}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Line-by-line modification of files}
\description{
Modifies a data file line by line, i.e. reads a file line by line, converts each line, then writes to the modified file. This method is especially useful when modifying large datasets, where the reading of entire files may be time consuming and require a large amount of memory.
}
\usage{
lineByLine(infile, outfile, linefunc = identity, choose.lines = NULL, 
choose.columns = NULL, ask = TRUE, blank.lines.skip = TRUE, verbose = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{infile}{A character string giving the name and path of the file to be modified. }
	\item{outfile}{A character string giving the name of the modified file. The name of the file is relative to the current working directory, unless the file name contains a definite path.}
	\item{linefunc}{\code{lineByLine} modifies each line using \code{linefunc}. Default is the identity function. The user may define his or her own line-modifying functions, see Details for a thorough description.}
	\item{choose.lines}{A numeric vector of lines to be selected from \code{infile}. If "NULL" (default), all lines are selected.}
	\item{choose.columns}{A numeric vector of columns to be selected from \code{infile}. By default, all columns are selected without reordering among the columns. Duplication and reordering among the selected columns will occur in the modified file corresponding to the order in which the columns are listed.}
	\item{ask}{Logical. Default is "TRUE". If set to "FALSE", an already existing outfile will be overwritten without asking.}
	\item{blank.lines.skip}{Logical. If "TRUE" (default), \code{lineByLine} ignores blank lines in the input.}
	\item{verbose}{Logical. Default is "TRUE", which means that the line number is displayed for each iteration, in addition to output from \code{linefunc}. If \code{choose.columns} contains invalid column numbers, this will also be displayed.}
	\item{...}{Further arguments to be passed to \code{linefunc}.}
}
\details{
	When reading large datafiles, functions such as \code{\link[utils]{read.table}} can use a large amount of memory and be extremely time consuming.
	Instead of reading the entire file at once, \code{lineByLine} reads one line at a time, modifies the line using \code{linefunc}, and then writes the line to \code{outfile}. \cr
		The user may specify his or her own line-converting function. This function must take the argument \code{x}, a character vector representing a single line of the file, split at spaces. However, additional arguments may be included.
	If \code{verbose} equals "TRUE", output should be displayed.
	The modified vector is returned. \cr	
	The framework of the line-modifying function should look something like this:
\preformatted{lineModify <- function(x){

.xnew <- x

## Define any modifications, for instance recoding missing values in 
## a dataset from NA to 0:
.xnew[is.na(.xnew)] <- 0

## Just to monitor progress, display, for instance, 10 first elements, 
## without newline:
cat(paste(.xnew[1:min(10, length(.xnew))], collapse = " "))

## Return converted vector
return(.xnew)
}}
	See \code{Haplin:::lineConvert} for an additional example of a line-modifying function.
}
	
\value{
  \code{lineByLine} returns the number of lines read, although invisible. The main objective is the modified file.
}

\references{
Web Site: \url{http://www.uib.no/smis/gjessing/genetics/software/haplin/}}

\author{Miriam Gjerdevik,\cr
with Hakon K. Gjessing\cr
Professor of Biostatistics\cr
Division of Epidemiology\cr
Norwegian Institute of Public Health\cr
\email{hakon.gjessing@fhi.no}}

\seealso{\code{\link[Haplin]{convertPed}}}

\examples{

\dontrun{

## Extract the first ten columns from "myfile.txt", without reordering
lineByLine(infile = "myfile.txt", outfile = "myfile_modified.txt", choose.columns = c(1:10))
}

}

