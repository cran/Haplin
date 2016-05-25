haptable.haplinStrat <- function(object){

## Create haptable for each stratum
.ut <- f.haptable.list(object)
## Change the generic name to stratum
names(.ut)[names(.ut) == "element"] <- "stratum"
## Change row.no to row.str to avoid confusion with row.win
names(.ut)[names(.ut) == "row.no"] <- "row.str"

return(.ut)
}
