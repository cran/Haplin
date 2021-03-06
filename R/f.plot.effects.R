f.plot.effects <- function(coeff, ref.cat, reference.method, haplos, markernames, maternal, poo, type = 1, ylim = c(0.2, 5), lwd = 2, use.single = seq(along = haplos), use.dd = seq(along = haplos), verbose = T, ...)
{
## PLOTS THE RESULT OF ESTIMATED ALLELE EFFECTS
##
#
## DECIDE WHETHER OR NOT TO USE THE ACTUAL HAPLOTYPE NAMES IN THE PLOT
.print.haplos <- T
## CHOOSE ROTATION OF HAPLOTYPE TEXT ON X-AXIS
.las <- 2 # 2: PERPENDICULAR TO AXIS, 0: ALWAYS PARALLEL TO AXIS
#
## INITIALIZE
.haplos <- haplos
.nall <- length(.haplos)
.mar <- par()$mar
.space <- 0.5
.xlab <- "Haplotype no."
if(.print.haplos) .xlab <- "Haplotype"
if(.print.haplos & (.las == 2)) .xlab <- "" # TURN OFF LABEL IF HAPLOTYPE NAMES ARE PERPEND. TO AXIS
.shift <- 0.05 # SEPARATION BETWEEN SINGLE AND DOUBLE DOSE IS 2*.shift
if(poo) .shift <- .shift * 2 ## NEED SPACE FOR TWO SINGLE DOSE EFFECTS
.len <- 0.02 # LENGHT OF CROSSBAR AT END OF CI IS 2*.len
.top <- T
.bottom <- T
#
## SELECT SINGLE- AND DOUBLE DOSES TO BE PLOTTED
.use.single <- use.single
.use.dd <- use.dd
if(reference.method == "ref.cat"){
	.use.single <- setdiff(.use.single, ref.cat)
	if(.nall == 2) .use.dd <- setdiff(.use.dd, ref.cat)
}


if(!missing(markernames)){
	if(length(markernames) == 1) .tmp <- "Marker: "
	else .tmp <- "Markers: "
	.markernames <- paste(.tmp, paste(markernames, collapse = "-"), sep = "")
	.width <- strwidth(.markernames, units = "figure")
}
#
### NB! FOLGENDE BOR VURDERES:
.yticks <- c(0.25, 0.5, 1, 2, 4)
#
##
if(type == 1){# CHILD ALONE
	.main <- "Relative risks for child haplotypes"
	.ylab <- "Relative risk (log scale)"
	.sel <- "c"
	if(.las == 2) .mar[1] <- .mar[1] + 1 # EXTEND LOWER MARGIN A LITTLE TO ACCOMMODATE LONG HAPLOTYPE NAMES
}
if(type == 2){# MOTHER ALONE
	.main <- "Relative risks for maternal haplotypes"
	.ylab <- "Relative risk (log scale)"
	.sel <- "m"
	if(.las == 2) .mar[1] <- .mar[1] + 1 # EXTEND LOWER MARGIN A LITTLE TO ACCOMMODATE LONG HAPLOTYPE NAMES
	if(!maternal) stop("Maternal effects must be estimated before they can be
    plotted!\n", call. = F)	#
}
if(type == 3){# CHILD, TOP HALF
	.main <- "Relative risks for haplotypes (log scale)"
	.ylab <- "Child"
	.sel <- "c"
	.mar[1] <- .space
	.top <- T
	.bottom <- F
}
if(type == 4){# MOTHER, BOTTOM HALF
	.main <- NULL
	.ylab <- "Mother"
	.sel <- "m"
	.mar[3] <- .space
	if(!maternal) stop("Maternal effects must be estimated before they can be
    plotted!\n", call. = F)	#
	.top <- F
	.bottom <- T
}
#
## PREPARE MESSAGE ABOUT REFERENCE
if(reference.method == "population") .ref.message <- "Ref = population"
else if (reference.method == "reciprocal") .ref.message <- "Ref = reciprocal"
else if(.print.haplos) .ref.message <- paste("Ref = ", names(ref.cat))
else .ref.message <- paste("Ref = ", ref.cat)
#
## SET MARGINS
par(mar = .mar)
#
## SET UP BASIC PLOT
if(.print.haplos){
	plot(1, 1, ..., xlim = c(0.5, .nall + 0.5), ylim = ylim, type = "n", xlab = .xlab, ylab = .ylab, log = "y", axes = F, main = .main, font = 2, font.lab = 2, xaxs = "i", yaxs = "i", cex.main = 1) #
	.haplos <- paste(.haplos, paste("(", round(100*coeff[1:.nall,1],1), "%)", sep = ""), sep = "\n")
	if(.bottom) axis(side = 1, at = 1:.nall, labels = .haplos, tick = F, font = 2, lwd = lwd, las = .las, cex.axis = 0.9)
}else{
	plot(1, 1, ..., xlim = c(0.5, .nall + 0.5), ylim = ylim, type = "n", xlab = .xlab, ylab = .ylab, log = "y", axes = F, main = .main, font = 2, font.lab = 2, cex.main = 1) #
	if(.bottom) axis(side = 1, at = 1:.nall, tick = F, font = 2, lwd = lwd)
}
#
## ADD MESSAGE ABOUT REFERENCE
if(!poo){
	.mtext <- '"s"'
}else{
	if(type == 1){# CHILD ALONE
		.mtext <- '"m" or "p"'
	}
	if(type == 2){# MOTHER ALONE
		.mtext <- '"s"'
	}
	if(type %in% 3:4){# CHILD TOP (OR MOTHER BOTTOM, NOT NEEDED)
		.mtext <- '"m", "p", or "s"'
	}
}
.mtext <- paste('Single dose = ', .mtext, '.  Double dose = "d".   ', .ref.message, '.', sep = "")
if(.top) mtext(.mtext, font = 2, cex = 0.8, line = 0.3)
#if(type == 3) .side <- 4#1
#if(type %in% 1:2) .side <- 4 
.at <- 1
if(type == 3) .at <- ylim[1]
if((type != 4) & !missing(markernames)) mtext(text = .markernames, side = 4, at = .at, cex = min(1, 0.7/.width), font = 2, col = grey(0.5))
#
## HORIZONTAL REFERENCE LINE
abline(h = 1, lwd = lwd)	#
#
## ADD Y-AXSIS
if(missing(ylim)) axis(side = 2, at = .yticks, tick = 0.02, font = 2, lwd = lwd)
else axis(side = 2, tick = 0.02, font = 2, lwd = lwd) 
#
## SETTING UP POSITIONS FOR BARS AND REFERENCE MARKER:
## SINGLE DOSE:
.pos <- 1:.nall - .shift #
.ref.pos <- .pos[ref.cat]
## SINGLE PATERNAL DOSE, IF POO:
.midpos <- 1:.nall
## DOUBLE DOSE:
.ddpos <- 1:.nall + .shift	#
#
## USE ONLY THOSE SELECTED
.pos <- .pos[.use.single]
.midpos <- .midpos[.use.single]
.ddpos <- .ddpos[.use.dd]
#
## MARK REFERENCE CATEGORY
if(reference.method == "ref.cat") textlabel(x = .ref.pos - .shift, y = 1, labels = "REF", cex = 0.7, font = 2)
### if(reference.method == "ref.cat") text(.ref.pos - .shift, 1 - 0.08, "REF", cex = 0.7, font = 2)	#
#
##
if(.sel == "m" | !poo){
	#
	## NAMES FOR RR PARAMETERS, ONLY USED TO EXTRACT VALUES FROM TABLE:
	.names <- paste("RR", .sel, 1:.nall, sep = "")
	.names <- .names[.use.single]
	#
	## EXTRACT COEFFICIENT VALUES
	.est <- coeff[.names, "est."]
	#
	.L <- coeff[.names, "lower"]
	.U <- coeff[.names, "upper"]
	#
	## BASIC PLOTTING OF SINGLE DOSE EFFECTS AND CIs
	.est.in <- f.Rplot(lwd = lwd, ylim = ylim, L = .L, U = .U, len = .len, pos = .pos, est = .est, pch = "s")
}
if(.sel == "c" & poo){
	#
	## NAMES FOR RR PARAMETERS, ONLY USED TO EXTRACT VALUES FROM TABLE:
	.names.mat <- paste("RR", "cm", 1:.nall, sep = "")
	.names.mat <- .names.mat[.use.single]
	#
	## EXTRACT COEFFICIENT VALUES
	.est.mat <- coeff[.names.mat, "est."]
	#
	.L.mat <- coeff[.names.mat, "lower"]
	.U.mat <- coeff[.names.mat, "upper"]
	#
	## BASIC PLOTTING OF SINGLE DOSE EFFECTS AND CIs
	.est.mat.in <- f.Rplot(lwd = lwd, ylim = ylim, L = .L.mat, U = .U.mat, len = .len, pos = .pos, est = .est.mat, pch = "m")
	#
	## NAMES FOR RR PARAMETERS, ONLY USED TO EXTRACT VALUES FROM TABLE:
	.names.pat <- paste("RR", "cf", 1:.nall, sep = "")
	.names.pat <- .names.pat[.use.single]
	#
	## EXTRACT COEFFICIENT VALUES
	.est.pat <- coeff[.names.pat, "est."]
	#
	.L.pat <- coeff[.names.pat, "lower"]
	.U.pat <- coeff[.names.pat, "upper"]
	#
	## BASIC PLOTTING OF SINGLE DOSE EFFECTS AND CIs
	.est.pat.in <- f.Rplot(lwd = lwd, ylim = ylim, L = .L.pat, U = .U.pat, len = .len, pos = .midpos, est = .est.pat, pch = "p")
	#
	## SAMLEVERDI, TIL RAPPORTERING
	.est.in <- .est.mat.in | .est.pat.in
}
#
## PLOT DOUBLE DOSE
.ddnames <- paste("RR", .sel, "dd", 1:.nall, sep = "")
.ddnames <- .ddnames[.use.dd]
#
.est.dd <- coeff[.ddnames, "est."]
#
.L.dd <- coeff[.ddnames, "lower"]
#
.U.dd <- coeff[.ddnames, "upper"]
#
## BASIC PLOTTING OF DOUBLE DOSE EFFECTS AND CIs
.est.dd.in <- f.Rplot(lwd = lwd, ylim = ylim, L = .L.dd, U = .U.dd, len = .len, pos = .ddpos, est = .est.dd, pch = "d")
#
## REMARK ABOUT PLOTTING RANGE
if(any(!.est.in) | any(!.est.dd.in)) {
	if(verbose) cat('\nNote: Some relative risk estimates fall outside the default plotting range.\nConsider replotting, with argument "ylim" set wider\n')
}
#
## RECTANGLE AROUND PLOT
rect(xleft = 0.5, ybottom = ylim[1], xright = .nall + 0.5, ytop = ylim[2], lwd = lwd)
#
##	
return(invisible(coeff))
}

