f.catch0 <- function(mcall, defaults){
## EXTRACT CALL
.mcall <- lapply(mcall[-1], function(x) eval.parent(x, 4))
#
## STANDARD haplin DEFAULTS
.defaults.hap <- formals(haplin0)
#
## OVERRIDE STANDARD haplin DEFAULTS
.defaults.hap[names(defaults)] <- defaults
#
## CHECK CALL FOR CONSISTENCY, SET OTHER DEFAULTS, ETC.
.info <- f.check.pars0(.mcall, .defaults.hap)
#####################
#
##
return(.info)
}
