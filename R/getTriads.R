
#' Getter for all full triads
#'
#' Wrapper function for \link{genDataGetPart} that returns a subset of the data containing 
#'  only full triads (where all, the child, the mother and the father have genetic data).
#'
#' @param data.in The data object (in format as the output of \link{genDataRead}); note
#'   that the design of the data is assumed to be triad.
#' @param file.out The base for the output filename (default: "my_data_onlyTriads").
#' @param dir.out The path to the directory where the output files will be saved.
#' @param overwrite Whether to overwrite the output files: if NULL (default), will prompt
#'   the user to give answer; set to TRUE, will automatically overwrite any existing files;
#'   and set to FALSE, will stop if the output files exist.
#'
#' @return A list object with three elements:
#'   \itemize{
#'     \item \emph{cov.data} - a \code{data.frame} with covariate data (if available in
#'        the input file)
#'     \item \emph{gen.data} - a list with chunks of the genetic data; the data is divided
#'        column-wise, using 10,000 columns per chunk; each element of this list is a
#'        \link[ff]{ff} matrix
#'     \item \emph{aux} - a list with meta-data and important parameters.
#'   }
#'   This now contains only the selected subset of data.
#'
#' @export
getFullTriads <- function(data.in = stop("No data given!", call. = FALSE),
	file.out = "my_data_onlyTriads",
	dir.out = ".", overwrite = NULL){

	if(!is(data.in, "haplin.data") ||
		!all(names(data.in) == .haplinEnv$.haplin.data.names)){
		stop("The input data is not in the correct format!", call. = FALSE)
	}
	
	format <- data.in$aux$info$filespecs$format
	
	if(format == "ped"){
		# first - get all the IDs and check family structure
		new.ids <- f.check.unique.ids(data.in$cov.data)
		id <- new.ids$ids
		pedIndex <- new.ids$pedIndex

  	# check that none of the children are missing
  	pedIndex <- new.ids$pedIndex <- pedIndex[!is.na(pedIndex[,"id.child"]), ]
  	# check if there are any families that have missing members in pedIndex
  	any.NAs.families <- apply(pedIndex, 1, function(row){
  	  any(is.na(row))
  	})
  	pedIndex <- new.ids$pedIndex <- pedIndex[!any.NAs.families, ]
	
		which.gen.data.fam <- f.create.missingness.matrix(data.in, new.ids)

		# here are IDs of the members where the full family information is available
		full.triads <- apply(which.gen.data.fam[ ,-1 ], 1, all)
		pedIndex.triads <- pedIndex[ full.triads, ]
		
		# check how many families found
		if(nrow(pedIndex.triads) == 0){
			stop("No full triads found!", call. = FALSE)
		}
		
		full.triads.rows <- sort(c(match(pedIndex.triads[ ,'id.father' ], id),
													match(pedIndex.triads[ ,'id.mother' ], id),
													match(pedIndex.triads[ ,'id.child' ], id)))
		return(genDataGetPart(data.in, design = "triad",
										rows = full.triads.rows,
										file.out = file.out, dir.out = dir.out, overwrite = overwrite)
					)
	} else if(format == "haplin"){
		stop("Not implemented yet", call. = FALSE)
	} else {
		stop(paste("Unrecognized format:", format), call. = FALSE)
	}
}

