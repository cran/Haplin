f.grid.asymp <- function(pos, design, xchrom, n.vars, nall, case.design, control.design){
#
	.design <- design
	.n.vars <- n.vars
	.case.design <- case.design
	.control.design <- control.design
	.pos <- pos
	#
	if(.design == "cc") .fam <- "c"
	else if(xchrom) .fam <- "mfx"
	else .fam <- "mf"
	if(.n.vars!=0) .pos <- .pos/(2*.n.vars)
	.grid <- f.pos.to.haplocomb(1:.pos, nall, fam = .fam)
	#
	## Add cc and sex columns
	if(xchrom) .grid <- rbind(cbind(sex=1,.grid),cbind(sex=2, .grid))
	if(.design!="triad") .grid <- rbind(cbind(cc=1,.grid),cbind(cc=2,.grid)) 
	#
	## Set columns missing (if family members are missing by design) and add children to .grid
	if(.design!="cc"){
		if(!xchrom){
			## Add children
			.grid.children <- lapply(1:length(nall), function(x){cbind(c1=.grid[,(length(nall)+x+ .n.vars)], c2=.grid[,(3*length(nall)+x+.n.vars)])})
			.grid.children <- do.call("cbind",.grid.children)
			colnames(.grid.children) <- paste("l",rep(1:length(nall),each=2),c(".c1",".c2"),sep="")
			.grid <- cbind(.grid,.grid.children)
		}		
		if(xchrom){
			## Add fathers and children
			.grid.fathers <- as.matrix(.grid[,grep("f",colnames(.grid))])
			colnames(.grid.fathers) <- paste("l",rep(1:length(nall),each=1),c(".f1"),sep="")
			.grid <- cbind(.grid,.grid.fathers)
			.grid.children <- lapply(1:length(nall), function(x){
				.c2 <- .grid[,(2*length(nall)+x+.n.vars)]
				.c2[which(.grid[,"sex"]==1)] <- .grid[which(.grid[,"sex"]==1),(length(nall)+x+.n.vars)]
				.grid.children <- cbind(c1=.grid[,(length(nall)+x+ .n.vars)], c2=.c2)
				.grid.children
			})
			.grid.children <- do.call("cbind",.grid.children)
			colnames(.grid.children) <- paste("l",rep(1:length(nall),each=2),c(".c1",".c2"),sep="")
			.grid <- cbind(.grid,.grid.children)
		}	
		if(.design=="cc.triad"){
			.grid.cases <- .grid[which(.grid[,"cc"]==2),]
			.grid.controls <- .grid[which(.grid[,"cc"]==1),]
		} else .grid.cases <- .grid
		#
		## Set columns missing (if family members are missing by design)
		if(!grepl("m", .case.design)) .grid.cases[,grep("m",colnames(.grid.cases))] <- NA 
		if(!grepl("f", .case.design)) .grid.cases[,grep("f",colnames(.grid.cases))] <- NA 
		if(.design=="cc.triad"){
			if(!grepl("m", .control.design)) .grid.controls[,grep("m",colnames(.grid.controls))] <- NA 
			if(!grepl("f", .control.design)) .grid.controls[,grep("f",colnames(.grid.controls))] <- NA 
			if(!grepl("c", .control.design)) .grid.controls[,which(colnames(.grid.controls)=="c")] <- NA
			.grid <- rbind(.grid.controls,.grid.cases)
		} else .grid <- .grid.cases
	}
	#
	## Grid order
		.grid.order <- paste("l",rep(1:length(nall),each=6),c(".m1", ".m2", ".f1", ".f2", ".c1",".c2"),sep="")
		if(.design=="cc") .grid.order <- paste("l",rep(1:length(nall),each=2),c(".c1",".c2"),sep="")
		if(.design=="triad" & xchrom) .grid <- .grid[,c("sex",.grid.order)]
		if(.design=="triad" & !xchrom) .grid <- .grid[,.grid.order]
		if(.design=="cc.triad" & !xchrom) .grid <- .grid[,c("cc",.grid.order)]
		if(.design=="cc.triad" & xchrom) .grid <- .grid[,c("cc","sex",.grid.order)]
		if(.design=="cc" & !xchrom) .grid <- .grid[,c("cc",.grid.order)]
		if(.design=="cc" & xchrom){ 
			.grid <- .grid[,c("cc","sex",.grid.order)]
			for(j in 1:length(nall)) .grid[which(.grid[,"sex"]==1),((2*j)+.n.vars)] <- .grid[which(.grid[,"sex"]==1),((2*j-1)+.n.vars)]
		}
	#	
	.grid <- as.matrix(.grid)
	mode(.grid) <- "character"
	#
	return(.grid)
	
}