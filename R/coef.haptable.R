coef.haptable <- function(tab){
##
## EXTRACT AND FORMAT A COEFFICIENT TABLE FROM A HAPTABLE.
## THE COEFFICIENT TABLE HAS THE SAME FORMAT AS THE RESULT OF USING
## summary.tri.glm(res)$effects ON A tri.glm OBJECT
## IN ADDITION SOME INFORMATION IS GIVEN AS ATTRIBUTES
##
#
## KEEP ONLY LINES CONTAINING HAPLOTYPES (USUALLY NOT NEEDED EXCEPT IN THE RARE CASES WHERE THERE ARE MORE MARKERS THAN HAPLOTYPES)
.tab <- tab[!is.na(tab$haplos),]
#
## DEDUCE VARIOUS INFORMATION
.maternal <- !is.null(.tab$RRm.est.)
.haplos <- .tab$haplos
.n.sel.haplo <- length(.haplos)
#
## EXTRACT REFERENCE INFORMATION FROM TABLE
.ref <- sort(unique(.tab$reference))## REMOVES ANY MISSING
if(identical(.ref, "reciprocal") | identical(.ref, "population")) {
###if((.ref[1] == "reciprocal") | (.ref[1] == "population")) {
	.reference.method <- .ref
	.ref.cat <- NA
} else
if(identical(.ref, c(" - ", "ref"))){
###if((length(.ref) == 2) && all(.ref == c(" - ", "ref"))){
	.reference.method <- "ref.cat"
	.ref.cat <- which(.tab$reference == "ref")
} else stop("Could not figure out reference. Something seems to be wrong with the haptable.")
#
## EXTRACT HAPLOTYPE FREQUENCIES
.p <- .tab[, c("haplofreq", "haplofreq.lower", "haplofreq.upper")]
.p <- cbind(.p, NA)
.colnavn <- c("est.", "lower", "upper", "p.value")
colnames(.p) <- .colnavn
rownames(.p) <- paste("p", 1:.n.sel.haplo, sep = "")
#
## EXTRACT CHILD RELATIVE RISK
.RR <- .tab[, c("RR.est.", "RR.lower", "RR.upper", "RR.p.value")]
colnames(.RR) <- .colnavn
rownames(.RR) <- paste("RRc", 1:.n.sel.haplo, sep = "")
#
.RRdd <- .tab[, c("RRdd.est.", "RRdd.lower", "RRdd.upper", "RRdd.p.value")]
colnames(.RRdd) <- .colnavn
rownames(.RRdd) <- paste("RRcdd", 1:.n.sel.haplo, sep = "")
#
## OUTPUT MATRIX
.ut <- rbind(.p, .RR, .RRdd)
#
## EXTRACT FOR MATERNAL, IF RELEVANT
if(.maternal){
	.RRm <- .tab[, c("RRm.est.", "RRm.lower", "RRm.upper", "RRm.p.value")]
	colnames(.RRm) <- .colnavn
	rownames(.RRm) <- paste("RRm", 1:.n.sel.haplo, sep = "")
	#
	.RRmdd <- .tab[, c("RRmdd.est.", "RRmdd.lower", "RRmdd.upper", "RRmdd.p.value")]
	colnames(.RRmdd) <- .colnavn
	rownames(.RRmdd) <- paste("RRmdd", 1:.n.sel.haplo, sep = "")
	#
	.ut <- rbind(.ut, .RRm, .RRmdd)
}
#
## ADD ATTRIBUTES
attr(.ut, "haplos") <- .haplos
attr(.ut, "maternal") <- .maternal
attr(.ut, "reference.method") <- .reference.method
attr(.ut, "ref.cat") <- .ref.cat
#
##
return(.ut)
}
