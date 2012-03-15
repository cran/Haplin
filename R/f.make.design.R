f.make.design <- function(maternal, response, info, test = F, ret.characteristics = F){
#
# THE PROGRAM ESTIMATES EFFECTS OF SEVERAL ALLELES IN A CASE-TRIAD, CASE-CONTROL-TRIAD OR CASE-CONTROL DESIGN
# CREATES AN APPROPRIATE DESIGN MATRIX FOR USE IN f.tri.glm
#
# THE MODEL ASSUMES HARDY-WEINBERG EQUILIBRIUM (AND RARE DISEASE, IF NECESSARY, IN CC)
# IMPORTANT: OBSERVED FREQUENCIES (observed) MUST COME FROM A COMPLETE
# GRID IN APPROPRIATE ORDER! THE ORDER IS:
# TRIAD: BY MOTHER'S FIRST ALLELE AND SECOND ALLELE,
# THEN THE FATHER'S FIRST AND SECOND
# CASE-CONTROL-TRIAD: BY MOTHER'S FIRST ALLELE AND SECOND ALLELE,
# THEN THE FATHER'S FIRST AND SECOND, AND FINALLY CONTROL (1) AND CASE (2), 
# SO THAT TOP HALF OF THE DATA ARE CONTROLS
# CASE-CONTROL: BY FIRST ALLELE, SECOND ALLELE
# AND FINALLY CONTROL (1) AND CASE (2), SO THAT TOP HALF OF THE DATA ARE
# CONTROLS
#
# maternal DETERMINES WHETHER OR NOT MATERNAL EFFECTS SHOULD BE ESTIMATED
# (NOT AVAILABLE FOR THE CASE-CONTROL DESIGN)
#
#
## EXTRACT VARIABLES
.ref.cat <- info$haplos$ref.cat
.xchrom <- info$model$xchrom
.design <- info$model$design
.sel.sex <- info$variables$sel.sex
.sel.sex.yes <- !is.null(.sel.sex)
.covar.yes <- !is.null(info$variables$covar)
.A <- sum(info$haplos$selected.haplotypes) ## NUMBER OF ALLELES (HAPLOTYPES)
#
###########
# CREATE "BASE" GRID OF PARENTAL GENOTYPES (AND POSSIBLY SEX AND/OR CASE-CONTROL STATUS)
###########
.tmp <- factor(1:.A)
.mf.tmp <- list(m1 = .tmp, m2 = .tmp, f1 = .tmp, f2 = .tmp, sex = c(0,1), cc = c(0,1))
## HERE, 0 AND 1 (MALES AND FEMALES, RESP.) ARE USED AS DUMMY FOR sex. ORIGINAL CODING (OUTSIDE f.tri.glm AND f.make.design) SHOULD BE 1 (MALES) AND 2 (FEMALES)
## HERE, 0 AND 1 (CONTROL AND CASE, RESP.) ARE USED AS DUMMY FOR cc. ORIGINAL CODING (OUTSIDE f.tri.glm AND f.make.design) SHOULD BE 1 (CONTROL) AND 2 (CASE). IN INPUT FILE, THE LARGEST IS ALWAYS THE CASE
.cond.triad <- (.design %in% c("triad", "cc.triad"))
.cond.f1 <- !.xchrom & (.design != "cc")
.cond.sex <- .xchrom & !.sel.sex.yes
.cond.cc <- (.design %in% c("cc", "cc.triad"))
#
## m1
if(!.cond.triad) .mf.tmp[["m1"]] <- NULL
## m2
## f1
if(!.cond.f1)  .mf.tmp[["f1"]] <- NULL
## f2
## sex
if(!.cond.sex) .mf.tmp[["sex"]] <- NULL
## cc
if(!.cond.cc) .mf.tmp[["cc"]] <- NULL
#
##
if((.design == "cc") & .xchrom) stop("Not implemented")## MERK: DE OVENSTÅENDE SELEKSJONENE SKAL FUNGERE FOR DENNE SITUASJONEN OGSÅ!



if(.covar.yes) .mf.tmp$covar <- factor(seq(along = info$variables$covar.codes))


if(ret.characteristics){
	.char <- sapply(.mf.tmp, length)
	## LITT AD-HOC, HÅPER DET IKKE TRENGS ANDRE STEDER:
	if(.design == "cc"){
		names(.char)[names(.char) == "m2"] <- "c1"
		names(.char)[names(.char) == "f2"] <- "c2"
	}
	return(.char)
}
#
## EXPAND GRID
.mf <- do.call("expand.grid", .mf.tmp)
#
###########
# CREATE DUMMY VARIABLES FOR PARENTAL GENOTYPES
###########
if(.cond.triad){
	.m1.dum <- model.matrix( ~ -1 + m1, data = .mf)
}
.m2.dum <- model.matrix( ~ -1 + m2, data = .mf)
if(.cond.f1){
	.f1.dum <- model.matrix( ~ -1 + f1, data = .mf)
}
.f2.dum <- model.matrix( ~ -1 + f2, data = .mf)
#
###########
# CREATE COUNTING VARIABLES FOR PARENTS AND CHILDREN
###########
## MOTHER
if(.cond.triad){
	.m.dumsum <- .m1.dum + .m2.dum
}else{
	.m.dumsum <- .m2.dum
}
## FATHER
if(.cond.f1){
	.f.dumsum <- .f1.dum + .f2.dum
}else{
	.f.dumsum <- .f2.dum
}
#
## PARENTS COMBINED, FOR ESTIMATING ALLELE FREQUENCIES
.parents.dumsum <- .m.dumsum + .f.dumsum
#
## CHILD, DEFINE .c.dumsum:
if(.cond.triad & .xchrom){
	## X-INACTIVATION(comb.sex == "double"): BOYS GET TWICE THE EFFECT. GIRLS GET ONE OR TWO.
	## MEANS THAT SINGLE DOSE IN BOYS CORRESPONDS TO DOUBLE DOSE IN GIRLS,
	## SINGLE DOSE IN GIRLS IS SQUARE ROOT OF DOUBLE DOSE (WHEN response = "mult")
	## BOYS: 1, RR^2, GIRLS: 1, RR, RR^2 (response = "mult")
	## BOYS: 1, RR1^2*RR2, GIRLS: 1, RR1, RR1^2*RR2 (response = "free")
	#
	## DOSE-RESPONSE(comb.sex == "single"): BOYS GET SINGLE DOSE EFFECT. GIRLS GET ONE OR TWO.
	## MEANS THAT SINGLE DOSE IN BOYS CORRESPONDS TO SINGLE DOSE IN GIRLS,
	## DOUBLE DOSE IN GIRLS IS SQUARE OF SINGLE DOSE. (WHEN response = "mult")
	## BOYS: 1, RR, GIRLS: 1, RR, RR^2 (response = "mult")
	## BOYS: 1, RR1, GIRLS: 1, RR1, RR1^2*RR2 (response = "free")
	#
	## DUMMIES FOR GIRLS
	.girls.dumsum <- .m2.dum + .f2.dum
	## DUMMIES FOR BOYS
	if(info$variables$comb.sex == "double"){
		## BOYS ARE COUNTED AS DOUBLE DOSE
		.boys.dumsum <- 2 * .m2.dum
		## WHEN BOTH BOYS AND GIRLS ARE INCLUDED, COMBINE
		.c.dumsum <- .boys.dumsum * (1- .mf$sex) + .girls.dumsum * .mf$sex
	}
	if(info$variables$comb.sex == "single"){
		## BOYS ARE COUNTED AS SINGLE DOSE
		.boys.dumsum <- .m2.dum
		## WHEN BOTH BOYS AND GIRLS ARE INCLUDED, COMBINE
		.c.dumsum <- .boys.dumsum * (1- .mf$sex) + .girls.dumsum * .mf$sex
	}
	if(info$variables$comb.sex == "males"){
		## BOYS ONLY. COUNTED AS SINGLE DOSE
		.boys.dumsum <- .m2.dum
		.c.dumsum <- .boys.dumsum
	}
	if(info$variables$comb.sex == "females"){
		## GIRLS ONLY
		.c.dumsum <- .girls.dumsum
	}
}else{
	## CHILD, STANDARD AUTOSOMAL MODEL
	.c.dumsum <- .m2.dum + .f2.dum
}
#
## SET APPROPRIATE COLUMN NAMES:
if(.cond.triad){
	dimnames(.m.dumsum)[[2]] <- paste("m", 1:.A, sep = "")
	dimnames(.f.dumsum)[[2]] <- paste("f", 1:.A, sep = "")
}
dimnames(.c.dumsum)[[2]] <- paste("c", 1:.A, sep = "")
dimnames(.parents.dumsum)[[2]] <- paste("mf", 1:.A, sep = "")
#
## LET THE FORMULA APPLY ONLY FOR CASES.
## (WARNING! RARE DISEASE ASSUMPTION HERE!)
if(.cond.cc){
	.c.dumsum <- .c.dumsum * .mf$cc
}
if(.design == "cc.triad"){
	.m.dumsum <- .m.dumsum * .mf$cc
	.f.dumsum <- .f.dumsum * .mf$cc
}

#warning("lklsdkjfoeifjoe")
#.mf$sex <- .mf$sex * .mf$cc
#.mf$sex <- 0



#
## THE FOLLOWING IS ONLY NECESSARY IF RESPONSE IS FREE
if(response == "free"){
	## SPECIAL CODING FOR HOMOZYGOTES:
	.c.dd <- (.c.dumsum == 2) + 0	# MODEL: R^2*Rstar
	if(.design == "cc.triad"){
		.c.dd <- .c.dd * .mf$cc
	}
	# SET APPROPRIATE COLUMN NAMES:
	dimnames(.c.dd)[[2]] <- paste("cdd", 1:.A, sep = "")	#
	#
	if(maternal){
		.m.dd <- (.m.dumsum == 2) + 0	# MODEL: R^2*Rstar FOR MAT. EFF
		if(.design == "cc.triad"){
			.m.dd <- .m.dd * .mf$cc			
		}
		# SET APPROPRIATE COLUMN NAMES:
		dimnames(.m.dd)[[2]] <- paste("mdd", 1:.A, sep = "")	
	}# END if(maternal)
}# END if(response == "free")



#
##
if(F){
## DETTE KAN JO KANSKJE KOMME INN SOM DOBBELTDOSE FOR "cc",
## HVIS ANT. HAPLOTYPER ER STØRRE ENN 2 (ELLER HVA SOM TRENGS?)
	.cc <- NULL # BARE TULL, MEN FOR AA UNNGAA FEILMELDING I R CMD check
	.mm.c.dobbeldose <- (.c.dumsum == 2) + 0	# MODEL: R^2*Rstar
	.mm.c.dobbeldose <- .mm.c.dobbeldose * .cc$cc
	dimnames(.mm.c.dobbeldose)[[2]] <- paste("cdd", 1:.A, sep = "")	
}











# SET UP FINAL DESIGN MATRIX, REMOVING REDUNDANT .ref.cat-ROW.
# NOTE: DIALLELIC SITUATION REQUIRES THE REMOVAL OF ONE DOUBLE DOSE PARAMETER:




.design.matrix <- cbind(sex = .mf$sex, cc = .mf$cc)


if(.covar.yes){
	.tmpmat <- cbind(.parents.dumsum, .mf[, c("covar"), drop = F])
	.form <- paste("mf", 1:.A, sep = "", collapse = " + ")
	.form <- paste("~ -1 + (", .form, "):covar", sep = "")
	.form <- formula(.form)
	.parents.dumsum <- model.matrix(.form, data = .tmpmat)
	cat("kontroller at sortering etc. blir rett!\n")
	.navn <- colnames(.parents.dumsum)
	.navn <- gsub(":", "_", .navn)

#	.navn <- gsub("_covar1", "", .navn)
#	cat("ad hoc!\n")

	colnames(.parents.dumsum) <- .navn

}



#
## ADD .parents.dumsum, FOR ESTIMATING HAPLOTYPE FREQUENCIES
.design.matrix <- cbind(.design.matrix, .parents.dumsum)
#
## IF response IS MORE THAN simple, ADD THE NECESSARY DESIGN VARIABLES



if(response != "simple"){
	.c.dumsum <- .c.dumsum[,  - .ref.cat, drop = F]
	.design.matrix <- cbind(.design.matrix, .c.dumsum)
	if(response == "free"){
		if(.A == 2) .c.dd <- .c.dd[,  - .ref.cat, drop = F]
		.design.matrix <- cbind(.design.matrix, .c.dd)
	}
	if(maternal){
		.m.dumsum <- .m.dumsum[,  - .ref.cat, drop = F]
		.design.matrix <- cbind(.design.matrix, .m.dumsum)
		if(response == "free"){
			if(.A == 2) .m.dd <- .m.dd[,  - .ref.cat, drop = F]
			.design.matrix <- cbind(.design.matrix, .m.dd)
		}
	}
}

.design.matrix <- as.data.frame(.design.matrix)

if(test & (response == "mult")){

	# DETTE ER EN EGEN TEST FOR AA PROEVE AA MATCHE RESULTATET FRA XLRT,
	# VED EN KJOERING
	# jj <- haplin(testdX1, .xchrom = T, sex = 1, n.vars = 2, markers = 1, verbose = T, printout = F, response = "mult", reference = 2, use.missing = T)
	# OG DERETTER SE PAA  exp(jj$result$result$coef)
	# MERK AT f.tri.glm KALLER f.make.design MED test = T

	.sex <- .design.matrix$sex
	.c1 <- .design.matrix$c1
	#
	.a1 <- (.c1 == 1) & (.sex == 0)
	.a2 <- (.c1 == 0) & (.sex == 1)
	.a3 <- (.c1 == 1) & (.sex == 1)
	.a4 <- (.c1 == 2) & (.sex == 1)
	#
	.design.matrix$a1 <- .a1 + 0
	.design.matrix$a2 <- .a2 + 0
	.design.matrix$a3 <- .a3 + 0
	.design.matrix$a4 <- .a4 + 0
	#
	.design.matrix$sex <- .design.matrix$c1 <- NULL
	#
	.msum <- .m1.dum[,2] + .m2.dum[,2]
	.fsum <- .f2.dum[,2]
	.mu <- paste(.msum, .fsum, sep = "-")
	.mu <- model.matrix(~ -1 + factor(.mu))
	colnames(.mu) <- paste("mu", seq(length.out = dim(.mu)[2]), sep = "")
	.design.matrix$mf1 <- .design.matrix$mf2 <- NULL
	.design.matrix <- cbind(.mu, .design.matrix)

}
#
##
return(.design.matrix)

}
