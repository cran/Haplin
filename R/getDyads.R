#' Getter only for all dyads (child and one parent)
#'
#' Wrapper function for \link{genDataGetPart} that returns a subset of the data containing 
#'  only dyads (where the child and only one parent have genetic data), i.e., not triads.
#'
#' @param data.in The data object (in format as the output of \link{genDataRead}); note
#'   that the design of the data is assumed to be "triad".
#' @param file.out The base for the output filename (default: "my_data_onlyDyads").
#' @param dir.out The path to the directory where the output files will be saved
#'   (default: ".", the current directory).
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
getDyads <- function(data.in = stop("No data given!", call. = FALSE),
	file.out = "my_data_onlyDyads",
	dir.out = ".",
	overwrite = NULL){

	if(!is(data.in, "haplin.data") ||
		!all(names(data.in) == .haplinEnv$.haplin.data.names)){
		stop("The input data is not in the correct format!", call. = FALSE)
	}
	
	format <- data.in$aux$info$filespecs$format
	
	if(format == "ped"){
		new.ids <- f.check.unique.ids(data.in$cov.data)
		id <- new.ids$ids
		pedIndex <- new.ids$pedIndex

  	# check that none of the children are missing
  	pedIndex <- new.ids$pedIndex <- pedIndex[!is.na(pedIndex[,"id.child"]), ]

  	# check which families have only one missing member in pedIndex
  	one.NA.families <- apply(pedIndex, 1, function(row){
  	  sum(is.na(row)) == 1
  	})

  	# if there is "NA" instead of a parental ID, this is a dyad
	  dyads.from.pedIndex <- c()
  	if (any(one.NA.families)){
    	dyads.from.pedIndex <- pedIndex[one.NA.families, ]
  	}
  	
  	# check if there are any families that have missing members in pedIndex
  	any.NAs.families <- apply(pedIndex, 1, function(row){
  	  any(is.na(row))
  	})
  	pedIndex <- new.ids$pedIndex <- pedIndex[!any.NAs.families, ]

  	# now - taking care of the genetic data missingness
  	which.gen.data.fam <- f.create.missingness.matrix(data.in, new.ids)

		# first, exclude the full triads
		full.triads <- apply(which.gen.data.fam[ ,-1 ], 1, all)
		pedIndex.no.triads <- pedIndex[ !full.triads, ]
		which.gen.data.fam <- which.gen.data.fam[ !full.triads, ]
		# then exclude the families with one member only
		single.mem.fam <- apply(which.gen.data.fam[ ,-1 ], 1,
			function(x){
				sum(x) == 1
				})
		pedIndex.dyads <- pedIndex.no.triads[ !single.mem.fam, ]
		which.gen.data.fam <- which.gen.data.fam[ !single.mem.fam, ]
		# remove the rows where only info from parents is available (if any such rows exist)
		only.parents <- apply(which.gen.data.fam[ ,-1 ], 1,
			function(x){
				! x[ "id.child" ]
				})
		if(any(only.parents)){
			pedIndex.dyads <- pedIndex.dyads[ !only.parents, ]
			which.gen.data.fam <- which.gen.data.fam[ !only.parents, ]
		}
		# check how many dyads found
		if (nrow(pedIndex.dyads) == 0 &
		    length(dyads.from.pedIndex) == 0){
			stop("No dyads found!", call. = FALSE)
		}
		if (nrow(pedIndex.dyads) == 0){
		  final.sel.IDs <- unique(
		    na.omit(as.character(
				  dyads.from.pedIndex[, -1]
				))
		 )
		} else if (length(dyads.from.pedIndex) == 0){
		  final.sel.IDs <- unique(c(
				pedIndex.dyads[ which.gen.data.fam$id.father ,'id.father' ],
				pedIndex.dyads[ which.gen.data.fam$id.mother ,'id.mother' ],
				pedIndex.dyads[ which.gen.data.fam$id.child ,'id.child' ]
			))
		} else {
		  final.sel.IDs <- unique(c(
				pedIndex.dyads[ which.gen.data.fam$id.father ,'id.father' ],
				pedIndex.dyads[ which.gen.data.fam$id.mother ,'id.mother' ],
				pedIndex.dyads[ which.gen.data.fam$id.child ,'id.child' ],
				na.omit(as.character(
				  dyads.from.pedIndex[, -1]
				))
			))
		}
		
		# find the correct rows, based on matching the IDs from the final selection
		only.dyads.rows <- sort(match(final.sel.IDs, id))
		return(genDataGetPart(data.in, design = "triad",
										rows = only.dyads.rows, file.out = file.out,
										dir.out = dir.out, overwrite = overwrite)
					)
	} else if(format == "haplin"){
		stop("Not implemented yet", call. = FALSE)
	} else {
		stop(paste("Unrecognized format:", format), call. = FALSE)
	}
}
