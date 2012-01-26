f.Rplot <- function(lwd, ylim, .L, .U, .L.dd, .U.dd, .len = .len, .pos, .est, .est.in, .ddpos, .est.dd, .est.dd.in, use.dd, use.single, ...){
##
## BASIC PLOTTING OF EFFECT VALUES AND THEIR CONFIDENCE INTERVALS
##
if(missing(use.dd)) use.dd <- seq(along = .est) # IF NOTHING ELSE REQUESTED, USE ALL


## LA DENNE INN FOR AA KORRIGERE FEIL I PLOTTING AV haptable. BURDE VAERT SKIKKELIG IMPLEMENTERT
if(missing(use.single)) use.single <- seq(along = .pos)


#
## CHECKING AND FIXING DATA POINTS OUTSIDE ylim. THIS REMOVES BOTHERSOME WARNINGS
## THAT DIDN'T DISAPPEAR WHEN SETTING suppress.graphics.warnings OR OTHER PARAMS.
#
#
.f.segments.rest <- function(x0, y0, x1, y1, y.lim, use = T, ...){
	## PLOTS ONLY STUFF WITHIN
	.above <- (y0 > y.lim[2]) & (y1 > y.lim[2])
	.below <- (y0 < y.lim[1]) & (y1 < y.lim[1])
	.use <- !.above & !.below & use
	
	
	y0 <- pmin(pmax(y.lim[1], y0), y.lim[2])
	y1 <- pmax(pmin(y.lim[2], y1), y.lim[1])
	
	if(sum(.use) == 0) return()
	segments(x0[.use], y0[.use], x1[.use], y1[.use], ...)
}


abline(h = 1, lwd = lwd)	#

#
## PLOTTING RELATIVE RISKS
points(.pos[.est.in], .est[.est.in], pch = "x", font = 2)	#
points(.ddpos[.est.dd.in], .est.dd[.est.dd.in], pch = "o", font = 2)	#
# WITH CIs AND LINE-ENDS
.f.segments.rest(.pos, .L, .pos, .U, y.lim = ylim, lwd = lwd, use = is.element(seq(along = .pos), use.single), ...)	#
.f.segments.rest(.pos - .len, .U, .pos + .len, .U, y.lim = ylim, lwd = lwd, use = is.element(seq(along = .pos), use.single), ...)
.f.segments.rest(.pos - .len, .L, .pos + .len, .L, y.lim = ylim, lwd = lwd, use = is.element(seq(along = .pos), use.single), ...)	#
#
.f.segments.rest(.ddpos, .L.dd, .ddpos, .U.dd, y.lim = ylim, lwd = lwd, use = is.element(seq(along = .ddpos), use.dd), ...)	#
.f.segments.rest(.ddpos - .len, .U.dd, .ddpos + .len, .U.dd, y.lim = ylim, lwd = lwd, use = is.element(seq(along = .ddpos), use.dd), ...)
.f.segments.rest(.ddpos - .len, .L.dd, .ddpos + .len, .L.dd, y.lim = ylim, lwd = lwd, use = is.element(seq(along = .ddpos), use.dd), ...)

}
