#' Checking uniqueness of individuals' IDs
#'
#' This checks whether the individuals have unique IDs. If not, it creates new IDs by
#'   adding the family ID to the indiv. ID. Since the function checks for specific
#'   column names, this basically means that the input should be as read from 
#'   a PED file.
#'
#' @param data.cov A matrix with "id.fam", "id.c", "id.m", and "id.f" columns.
#'
#' @return List with:
#'   \itemize{
#'     \item ids - new IDs in the same order as "id.c" column in the input data
#'     \item pedIndex - matrix with "family", "id.child", "id.mother", and
#'        "id.father" columns
#'   }
#'
#' @keywords internal

f.check.unique.ids <- function( data.cov ){
	# Basic checks. NOTE: id.c does not only refer to the child, but to all individuals in the file
	if(!identical(colnames(data.cov)[1:4], c("id.fam", "id.c", "id.f", "id.m"))) stop("Something's wrong with the id variables", call. = F)
	if(any(is.na(data.cov[,"id.fam"]) | is.na(data.cov[, "id.c"]))) stop("Missing individual id and/or family id", call. = F)
	#
	id.c <- data.cov[ ,"id.c" ]
	if( anyDuplicated(id.c) > 0 ){ # dupl mother, father is presumably OK
		cat( "   Creating unique IDs for individuals...\n" )
		tmp.data.cov <- data.cov[, 1:4]
		# combine fam and c id's
		tmp.data.cov[, "id.c"] <- paste(data.cov[, "id.fam"], data.cov[, "id.c"], sep = "_")
		# combine fam and f id's, but keep missing as missing
		tmp.data.cov[, "id.f"] <- paste(data.cov[, "id.fam"], data.cov[, "id.f"], sep = "_")
		tmp.data.cov[, "id.f"][is.na(data.cov[, "id.f"])] <- NA
		# combine fam and m id's, but keep missing as missing
		tmp.data.cov[, "id.m"] <- paste(data.cov[, "id.fam"], data.cov[, "id.m"], sep = "_")
		tmp.data.cov[, "id.m"][is.na(data.cov[, "id.m"])] <- NA
		#
		data.cov[, 1:4] <- tmp.data.cov
		#
		cat( "   ...done.\n" )
	} # end if any duplicated
	#
	id.c.new <- data.cov[ ,"id.c" ]
	if( anyDuplicated(id.c.new) > 0 ){
		#
		id.c.new.dupl <- id.c.new[duplicated(id.c.new)]
		print(id.c.new.dupl[1:min(4, length(id.c.new.dupl))], "...")
		stop("Found duplicated children's id's within same family", call. = F)
	}
	#
	# sort the families and check coding
	pedIndex <- f.prep.pedIndex( data.cov )
	#
	return( list( ids = id.c.new, pedIndex = pedIndex ) )
}
