f.insert.vector <- function(target, insert, pos.target, len.insert){
#
## INSERTS PIECES OF insert INTO target. THE PIECES ARE INSERTED STARTING
# IN POSITION pos.target. THE VALUES IN target ARE SHIFTED TO THE RIGHT.
# THE PIECES OF insert ARE OF LENGTH len.insert
# pos.target AND len.insert MUST BE OF THE SAME LENGTH.
# TO PLACE A PIECE AT THE END, USE A POSITION EQUAL TO length(target) + 1
#
	.n.target <- length(target)
	.n.insert <- length(insert)
	if(length(pos.target) != length(len.insert)) stop("Problem in vector insertion...")
	#
	#
	f.vis(.len.target <- diff(c(1, pos.target, .n.target + 1)), vis = F)
	f.vis(.len.insert <- c(len.insert, 0), vis = F)
	if(.n.target != sum(.len.target) | .n.insert != sum(.len.insert)) stop("Problem in vector insertion...")
	if(any(.len.target < 0) | any(.len.insert < 0)) stop("Problem in vector insertion...")
	
	.n.blocks <- length(.len.target)

	.ut <- c(target, insert) # CORRECT LENGTH AND MODE
	.ut[] <- NA
#
#
	.len <- as.numeric(t(matrix(c(.len.target, .len.insert), ncol = 2)))
	.ind <- rep(rep(1:2, .n.blocks), .len)

	.ut[.ind == 1] <- target
	.ut[.ind == 2] <- insert
	.ut
	}