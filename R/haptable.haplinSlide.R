haptable.haplinSlide <- function(object){

## Create haptable for each window
.ut <- f.haptable.list(object)
## change the generic name to window
names(.ut)[names(.ut) == "element"] <- "window"
## Change row.no to row.win to avoid confusion with row.str
names(.ut)[names(.ut) == "row.no"] <- "row.win"

return(.ut)
}
