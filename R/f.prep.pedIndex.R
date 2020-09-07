f.prep.pedIndex <- function( data.cov ){
## EXTRACT PED-INFORMATION

## EXTRACT FAMILY INFORMATION FROM PED FILE
# cat("Extracting family, sex and case/control information from ped file...\n")
.fam <- as.dframe( data.cov )

names( .fam ) <- c( "family", "id", "father", "mother", "sex", "cc" )
.pheno <- as.matrix( .fam[ , c( "id", "sex", "cc" ) ], mode = "character" )
.fam <- .fam[, c( "family", "id", "father", "mother" ) ]

.sex.u <- sort( unique( .pheno[ ,"sex" ] ) )
if( length( .sex.u ) > 2 ){
	stop( "More than 2 different codes in the sex column (column 5)", call. = F )
}
if( any( !is.element( .sex.u, c("0", "1") ) ) ){
	if( all( is.element( .sex.u, c("1", "2") ) ) ){
		.pheno[ .pheno[ ,"sex" ] == "2","sex" ] <- "0" # RECODE FEMALES
	}else{
		stop("Invalid codes in the sex column (column 5)", call. = F)
	}
}
.pheno <- as.dframe( .pheno )

## CREATE INDEXING TO BE USED LATER WHEN CONVERTING TO HAPLIN FORMAT
.pedIndex <- f.make.index( .fam, output = "ids" )

# return( invisible() )
return( .pedIndex )
}
