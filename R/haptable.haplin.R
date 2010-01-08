haptable.haplin <- function(object){
## CREATES A SUMMARY TABLE OF ESTIMATION RESULTS FROM HAPLIN
## object IS THE RESULT FROM A HAPLIN RUN
##
#
## COMPUTE SUMMARY
.summ.res <- summary(object)
#
## NUMBER OF TRIADS USED IN ANALYSIS
.ntri.seq <- object$ntri.seq
#.ntri.used <- as.numeric(object$ntri.seq["After rem rare haplos"])
#if(.ntri.used != object$result$ntri) stop() ## object$result$ntri ER IKKE ALLTID *HELT* AVRUNDET
#
## EXTRACT AND COMPUTE ALLELE-RELATED RESULTS
.alls <- t(sapply(object$alleles, function(x) c(alleles = paste(names(x), collapse = "/"), counts = paste(x, collapse = "/"))))
.alls <- data.frame(.alls, HWE.pv = sapply(object$HWE.res, function(x) x$p.value), stringsAsFactors = F)
#
## EXTRACT HAPLOTYPES USED
.selected.haplotypes <- object$selected.haplotypes
.nh <- sum(.selected.haplotypes)
.selected.haplotypes <- names(.selected.haplotypes)[.selected.haplotypes]
#
## EXTRACT EFFECTS MATRIX
.effs <- .summ.res$summary.tri.glm$effects
#
## REFORMAT TABLE
.tab <- data.frame(haplofreq = .effs[1:.nh,], RR = .effs[(.nh+1):(2*.nh),], RRdd = .effs[(2*.nh+1):(3*.nh),])
if(object$result$maternal){
	.tabm <- data.frame(RRm = .effs[(3*.nh+1):(4*.nh),], RRmdd = .effs[(4*.nh+1):(5*.nh),])
	.tab <- cbind(.tab, .tabm)
}
#
## JOIN RESULTS INTO DATA FRAME
.tab <- data.frame(haplos = .selected.haplotypes, pv.overall = rep(object$loglike["p.value.overall"], dim(.tab)[1]), .tab, stringsAsFactors = F)
.tab$haplofreq.p.value <- NULL


.diff <- .nh - dim(.alls)[1]
if(.diff > 0) .alls <- as.data.frame(lapply(.alls, function(x) c(x, rep(NA, .diff))))
if(.diff < 0) .tab <- as.data.frame(lapply(.tab, function(x) c(x, rep(NA, -.diff))))

#
##
.ntri.mat <- matrix(.ntri.seq, nrow = dim(.tab)[1], ncol = 4, byrow = T, dimnames = list(NULL, names(.ntri.seq)))

.tab <- data.frame(.ntri.mat, .alls, .tab)

names(.tab)[names(.tab) == "haplofreq.est."] <- "haplofreq"
rownames(.tab) <- NULL
return(.tab)

write.table(.tab, file = "test.dat", sep = "\t", quote = F, row.names = F) 

}

