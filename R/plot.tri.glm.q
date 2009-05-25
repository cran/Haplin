"plot.tri.glm"<- function(res, reference.method, design, selected.haplotypes, type = 1, ylim = c(0.2, 5), lwd = 2, use.dd, info, ...)
{
# PLOTS THE RESULT OF ESTIMATED ALLELE EFFECTS
#

.mar <- par()$mar
.space <- 0.5
#
## DECIDE WHETHER OR NOT TO USE THE ACTUAL HAPLOTYPE NAMES IN THE PLOT
.print.haplos <- T
## CHOOSE ROTATION OF HAPLOTYPES ON X-AXIS
.las <- 2 # 2: PERPENDICULAR TO AXIS, 0: ALWAYS PARALLEL TO AXIS
#
## INITIALIZE
#
.nall <- res$nall
.ref.cat <- res$ref.cat #
.haplos <- names(selected.haplotypes)[selected.haplotypes]
if(missing(use.dd)) use.dd <- seq(along = .haplos)
.maternal <- res$maternal
.xlab <- "Haplotype no."
if(.print.haplos) .xlab <- "Haplotype"
if(.print.haplos & (.las == 2)) .xlab <- "" # TURN OF LABEL IF HAPLOTYPE NAMES ARE PERPEND. TO AXIS
.shift <- 0.05
.len <- 0.02	#
.top <- T
.bottom <- T
#
### NB! FOLGENDE BOR VURDERES:
.yticks <- c(0.25, 0.5, 1, 2, 4)
#
#
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
	if(!.maternal) stop("Maternal effects must be estimated before they can be plotted!\n")	#
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
	if(!.maternal) stop("Maternal effects must be estimated before they can be plotted!\n")	#
	.top <- F
	.bottom <- T
}
#
## SET MARGINS
par(mar = .mar)
#
#
## EXTRACT VALUES FOR PLOTTING
#
.coef <- summary(res, design = design, reference.method = reference.method, conf.int = T, info = info)$effects #
#
#
## MARKING REFERENCE AND SETTING UP POSITIONS FOR BARS:
	.pos <- 1:.nall - .shift #
	.ref.pos <- .pos[.ref.cat]
	if(reference.method == "ref.cat").pos <- .pos[-.ref.cat]
	.ddpos <- 1:.nall + .shift	#
	if(reference.method == "ref.cat" & .nall == 2) .ddpos <- .ddpos[-.ref.cat]
#
## NAMES FOR RR PARAMETERS, ONLY USED TO EXTRACT VALUES FROM TABLE:
	.names <- paste("RR", .sel, 1:.nall, sep = "")
	if(reference.method == "ref.cat") .names <- .names[ - .ref.cat]
	.ddnames <- paste("RR", .sel, "dd", 1:.nall, sep = "")
	if(reference.method == "ref.cat" & .nall == 2) .ddnames <- .ddnames[ - .ref.cat] #
#
#
## SET UP BASIC PLOT:

if(reference.method == "population") .ref.message <- "Ref = population"
else if (reference.method == "reciprocal") .ref.message <- "Ref = reciprocal"
else if(.print.haplos) .ref.message <- paste("Ref = ", .haplos[.ref.cat])
else .ref.message <- paste("Ref = ", .ref.cat)
#
.f.in <- function(x, yl) {(x >= yl[1]) & (x <= yl[2])}
#
.est <- .coef[.names, "est."]
.est.in <- .f.in(.est, ylim)
.est.dd <- .coef[.ddnames, "est."]
.est.dd.in <- .f.in(.est.dd, ylim) & is.element(seq(along = .est.dd), use.dd) # PLOT ONLY EFFECTS WITHIN BOUNDARIES, AND WHICH THE USER REQUESTS
#
.L <- .coef[.names, "lower"]
.L.dd <- .coef[.ddnames, "lower"]
#
.U <- .coef[.names, "upper"]
.U.dd <- .coef[.ddnames, "upper"]
#

if(.print.haplos){
		plot(1, 1, ..., xlim = c(0.5, .nall + 0.5), ylim = ylim, type = "n", xlab = .xlab, ylab = .ylab, log = "y", axes = F, main = .main, font = 2, font.lab = 2, xaxs = "i", yaxs = "i", cex.main = 1) #
		#.underline <- c("\n^-^-^", "", "", "")
		#.a1 <- c("\n", "", "\n", "")
		#.a2 <- c("", "\n", "", "\n")
		#.haplos <- paste(.a1, .haplos, .haplos, .haplos, .a2, sep = "")
		
###		tull <<- .haplos
		.haplos <- paste(.haplos, paste("(", round(100*.coef[1:.nall,1],1), "%)", sep = ""), sep = "\n")
		if(.bottom) axis(side = 1, at = 1:.nall, labels = .haplos, tick = F, font = 2, lwd = lwd, las = .las, cex.axis = 0.9)
	}
	else{
		plot(1, 1, ..., xlim = c(0.5, .nall + 0.5), ylim = ylim, type = "n", xlab = .xlab, ylab = .ylab, log = "y", axes = F, main = .main, font = 2, font.lab = 2, cex.main = 1) #
		if(.bottom) axis(side = 1, at = 1:.nall, tick = F, font = 2, lwd = lwd)
	}


f.vis(.est.in, vis = F)
f.vis(.est.dd.in, vis = F)
if(any(!.est.in) | any(!.est.dd.in)) {
	if(info$control$verbose) cat('\nNote: Some relative risk estimates fall outside the default plotting range.\nConsider replotting, with argument "ylim" set wider\n')
}

f.Rplot(lwd = lwd, ylim = ylim, .L = .L, .U = .U, .L.dd = .L.dd, .U.dd = .U.dd, .len = .len, .pos = .pos, .est = .est, .est.in = .est.in, .ddpos = .ddpos, .est.dd = .est.dd, .est.dd.in = .est.dd.in, use.dd = use.dd)

	if(missing(ylim)) axis(side = 2, at = .yticks, tick = 0.02, font = 2, lwd = lwd)
	else
	axis(side = 2, tick = 0.02, font = 2, lwd = lwd) 

.mtext <- paste("Single dose = \"x\", Double dose = \"o\", ", .ref.message)

if(.top) mtext(.mtext, font = 2, cex = 0.8, line = 0.3)
if(reference.method == "ref.cat")text(.ref.pos - .shift, 1 - 0.08, "REF", cex = 0.7, font = 2)	#

if(F){
abline(h = ylim, lwd = 2)
segments(x0 = c(0, .nall + 0.5), x1 = c(0, .nall + 0.5), y0 = c(0,0), y1 = ylim, lwd = 2)
box(which = "outer")
}
rect(xleft = 0.5, ybottom = ylim[1], xright = .nall + 0.5, ytop = ylim[2])

#
#
#	
	return(invisible(.coef))

}
