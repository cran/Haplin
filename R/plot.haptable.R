plot.haptable <- function(tab, separate.plots = F, filename, filetype = "png", use.dd, verbose = T, ...)
{
##
## PLOT A HAPTABLE
## MERK: DENNE HAR MYE FELLES MED plot.haplin, BURDE KANSKJE VAERT SAMKJOERTE
##
#
.coef <- coef.haptable(tab)
.n.sel.haplo <- length(attr(.coef, "haplos"))
.maternal <- attr(.coef, "maternal")
.ref.cat <- attr(.coef, "ref.cat")
.reference.method <- attr(.coef, "reference.method")
.haplos <- attr(.coef, "haplos")
#
if(missing(use.dd)) use.dd <- 1:.n.sel.haplo
#
## PRODUCE JPEG OR PNG, IF REQUESTED
.prod.file <- !missing(filename)
#
##
.params <- list(coeff = .coef, ref.cat = .ref.cat, reference.method = .reference.method, haplos = .haplos, maternal = .maternal, use.dd = use.dd, verbose = verbose, ...)

##
## MERK! HAR IKKE IMPLEMENTERT use.single, SOM BRUKES TIL Å PLOTTE BOYS ONLY 
## I X-CHROM. DEN ER HELLER IKKE IMPLEMENTERT I f.plot.effects. SE IMPLEMENTERING I 
## plot.haplin



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
		.params$type <- 4
		invisible(do.call("f.plot.effects", .params))
	}
	if(.maternal & separate.plots){
		.params$type <- 1
		invisible(do.call("f.plot.effects", .params))
		.params$type <- 2
		invisible(do.call("f.plot.effects", .params))
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
		.params$type <- 4
		invisible(do.call("f.plot.effects", .params))
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
		.params$type <- 2
		invisible(do.call("f.plot.effects", .params))
		dev.off()	
	}
}# END .prod.file


}

