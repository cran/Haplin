f.Rplot <- function(lwd, ylim, L, U, len, pos, est, pch, ...){
##
## BASIC PLOTTING OF EFFECT VALUES AND THEIR CONFIDENCE INTERVALS
##
#
## CHECKING AND FIXING DATA POINTS OUTSIDE ylim. THIS REMOVES BOTHERSOME WARNINGS
## THAT DIDN'T DISAPPEAR WHEN SETTING suppress.graphics.warnings OR OTHER PARAMS.
#
.f.segments.rest <- function(x0, y0, x1, y1, ...){
	## PLOTS ONLY STUFF WITHIN
	.above <- (y0 > ylim[2]) & (y1 > ylim[2])
	.below <- (y0 < ylim[1]) & (y1 < ylim[1])
	# SELECT THOSE WITH _SOMETHING_ WITHIN
	.use <- !.above & !.below
	# RESTRICT TO WITHIN ylim	
	y0 <- pmin(pmax(ylim[1], y0), ylim[2])
	y1 <- pmax(pmin(ylim[2], y1), ylim[1])
	# IF NOTHING TO PLOT
	if(sum(.use) == 0) return()
	# PLOTTING, RESTRICTED
	segments(x0[.use], y0[.use], x1[.use], y1[.use], lwd = lwd, ...)
}
#
## PLOTTING RELATIVE RISKS
# WITH CIs AND LINE-ENDS
## VERTICAL BAR
.f.segments.rest(pos, L, pos, U, ...)
## HORIZONTAL TOP AND BOTTOM BARS
.f.segments.rest(pos - len, U, pos + len, U, ...)
.f.segments.rest(pos - len, L, pos + len, L, ...)
#
## PLOT ONLY THOSE WITHIN
.f.in <- function(x, yl) {(x >= yl[1]) & (x <= yl[2])}
.est.in <- (est >= ylim[1]) & (est <= ylim[2])
#
## EFFECT
points(pos[.est.in], est[.est.in], pch = 22, cex = 2, bg = "white", col = "white") # PROVIDES WHITE BACKGROUND FOR PLOTTING CHARACTER
points(pos[.est.in], est[.est.in], pch = pch, font = 2) # PLOTTING CHARACTER, SUCH AS "s", "d" ETC.
#
## RETURN WHICH ONES WERE OUTSIDE RANGE---- TO REPORT/RECOMMEND REPLOTTING
return(.est.in)
}
