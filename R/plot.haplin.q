plot.haplin <- function(x, reference, separate.plots = F, filename, filetype = "png", use.dd, ...)
{
#
.n.sel.haplo <- sum(x$selected.haplotypes)
.maternal <- x$result$maternal
.info <- x$info
.sel.sex <- .info$variables$sel.sex
.response <- .info$haplos$response
.poo <- .info$model$poo
.haplos <- names(x$selected.haplotypes)[x$selected.haplotypes]
.verbose <- .info$control$verbose
if(.info$model$scoretest == "only") stop('Sorry, plotting of effect estimates not available when haplin was run with scoretest = "only"')
if(!is.null(.sel.sex) && (.sel.sex == 1)){
	## IF SELECT ONLY BOYS, SHOW ONLY SINGLE DOSE
	.use.single <- 1:.n.sel.haplo
	.use.dd <- numeric(0)
}else{
	if(missing(use.dd)) .use.dd <- 1:.n.sel.haplo
	else .use.dd <- use.dd
	.use.single <- 1:.n.sel.haplo
}
if(missing(use.dd)) .use.dd.mat <- 1:.n.sel.haplo
else .use.dd.mat <- use.dd
.use.single.mat <- 1:.n.sel.haplo
#
## USE APPROPRIATE REFERENCE
if(missing(reference)){
	.reference.method <- x$reference.method
} 
else if (is.numeric(reference)){
	cat("\nWARNING: REFERENCE CATEGORY CAN ONLY BE SET IN FIRST RUN OF HAPLIN!\nFOR summary AND plot METHODS ONLY REFERENCE METHOD CAN BE CHOSEN, NOT CATEGORY\n\n")
	.reference.method <- x$reference.method
}
else .reference.method <- reference
#
## CHECK THAT ONLY REFCAT IS USED WHEN ONLY TWO HAPLOTYPES/ALLELES	
if(.n.sel.haplo == 2 & .reference.method != "ref.cat"){
	cat("\nNOTE: ONLY SINGLE REFERENCE CATEGORY METHOD ALLOWED FOR TWO HAPLOTYPES/ALLELES!\n (reference has been set to", x$result$ref.cat, ")\n")
	.reference.method <- "ref.cat"
}
#
## PRODUCE JPEG OR PNG, IF REQUESTED
.prod.file <- !missing(filename)
#
## EXTRACT COEFFICIENTS
.coef <- summary(x$result, reference.method = .reference.method, conf.int = T, info = .info)$effects
.ref.cat <- x$result$ref.cat
#
## PARAMETERS FOR CHILD PLOT
.params <- list(coeff = .coef, ref.cat = .ref.cat, reference.method = .reference.method, haplos = .haplos, maternal = .maternal, poo = .poo, use.dd = .use.dd, use.single = .use.single, verbose = .verbose, ...)
#
## CHANGE SOME PARAMETERS FOR MATENAL PLOT
.params.mat <- .params
.params.mat$use.dd <- .use.dd.mat
.params.mat$use.single <- .use.single.mat
#
##
if(!.prod.file){
	## RETAIN OLD PARAMETERS
	.oldpar <- par(no.readonly = T)
	on.exit(par(.oldpar))
	#
	if(!.maternal){
		.params$type <- 1
		invisible(do.call("f.plot.effects", .params))
	}
	if(.maternal & !separate.plots){
		par(mfrow = c(2,1), oma = c(2,0,0,0))
		.par <- par(no.readonly = T)
		.params$type <- 3
		invisible(do.call("f.plot.effects", .params))
		par(mar = .par$mar)
		.params.mat$type <- 4
		invisible(do.call("f.plot.effects", .params.mat))
	}
	if(.maternal & separate.plots){
		.params$type <- 1
		invisible(do.call("f.plot.effects", .params))
		.params.mat$type <- 2
		invisible(do.call("f.plot.effects", .params.mat))
	}
}# END !.prod.file


if(.prod.file){
	.jpeg.size <- c(440, 460)
	if(.n.sel.haplo > 4) .jpeg.size <- .jpeg.size + (.n.sel.haplo - 4) * c(30,0)
	#
	if(!.maternal){
		if(filetype == "png"){
			png(filename = filename, width = .jpeg.size[1], height = .jpeg.size[2], pointsize = 9)
		}
		if(filetype == "jpeg"){
			jpeg(filename = filename, width = .jpeg.size[1], height = .jpeg.size[2], pointsize = 9, quality = 100)
		}
		.params$type <- 1
		invisible(do.call("f.plot.effects", .params))
		dev.off()
	}
	if(.maternal & !separate.plots){
		if(filetype == "png"){
			png(filename = filename, width = .jpeg.size[1], height = .jpeg.size[2], pointsize = 9)
		}
		if(filetype == "jpeg"){
			jpeg(filename = filename, width = .jpeg.size[1], height = .jpeg.size[2], pointsize = 9, quality = 100)
		}
		par(mfrow = c(2,1), oma = c(2,0,0,0))
		.par <- par(no.readonly = T)
		.params$type <- 3
		invisible(do.call("f.plot.effects", .params))
		par(mar = .par$mar)
		.params.mat$type <- 4
		invisible(do.call("f.plot.effects", .params.mat))
		dev.off()
	}
	if(.maternal & separate.plots){
		if(filetype == "png"){
			png(filename = paste("1", filename, sep = ""), width = .jpeg.size[1], height = .jpeg.size[2], pointsize = 9)
		}
		if(filetype == "jpeg"){
			jpeg(filename = paste("1", filename, sep = ""), width = .jpeg.size[1], height = .jpeg.size[2], pointsize = 9, quality = 100)
		}
		.params$type <- 1
		invisible(do.call("f.plot.effects", .params))
		dev.off()
		if(filetype == "png"){
			png(filename = paste("2", filename, sep = ""), width = .jpeg.size[1], height = .jpeg.size[2], pointsize = 9)
		}
		if(filetype == "jpeg"){
			jpeg(filename = paste("2", filename, sep = ""), width = .jpeg.size[1], height = .jpeg.size[2], pointsize = 9, quality = 100)
		}
		.params.mat$type <- 2
		invisible(do.call("f.plot.effects", .params.mat))
		dev.off()	
	}
}# END .prod.file
}
