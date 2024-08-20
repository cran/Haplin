context("Testing preprocessing of data")

source("common_vars.R")
gen.data.read <- suppressWarnings(
  genDataRead(
    file.in = my.ped.data, format = "ped", file.out = my.ped.out.data,
    overwrite = if.overwrite
 )
)
gen.data.read <- suppressWarnings(
  genDataRead(
    file.in = my.haplin.data.no.cov, format = "haplin", n.vars = 0,
    file.out = my.haplin.no.cov.out, overwrite = if.overwrite
  )
)
gen.data.read <- suppressWarnings(
  genDataRead(
    file.in = my.haplin.data.cov, format = "haplin", n.vars = 2, allele.sep = "",
    col.sep = "\t", file.out = my.haplin.cov.out, overwrite = if.overwrite
 )
)
rm(gen.data.read)

test_that("Preprocessing haplin data", {
	haplin.data.read <- genDataLoad(my.haplin.cov.out)
	haplin.data.prep <- suppressWarnings(
	  genDataPreprocess(
	    data.in = haplin.data.read, file.out = my.haplin.cov.preproc.out,
	    overwrite = if.overwrite
	  )
	)
	
	expect_s3_class(haplin.data.prep, "haplin.ready")
	expect_equal(ncol(haplin.data.prep$gen.data[[1]]), ncol(haplin.data.read$gen.data[[1]]))
	expect_equal(nrow(haplin.data.prep$gen.data[[1]]), nrow(haplin.data.read$gen.data[[1]]))
})

test_that("Preprocessing ped data", {
	ped.data.read <- genDataLoad(my.ped.out.data)
	# I need to take only some markers, otherwise it would run for too long!
	ped.data.part <- suppressWarnings(
	  genDataGetPart(
	    ped.data.read, design = "triad", markers = 110:111, overwrite = if.overwrite
	  )
	)
	ped.data.prep <- suppressWarnings(
	  genDataPreprocess(
	    data.in = ped.data.part, file.out = my.ped.preproc.out, overwrite = if.overwrite
	  )
	)
	
	expect_s3_class(ped.data.prep, "haplin.ready")
	expect_equal(ncol(ped.data.prep$gen.data[[1]]), ncol(ped.data.part$gen.data[[1]]) * 3)
	expect_equal(nrow(ped.data.prep$gen.data[[1]]), 559)
	
	expect_equal(unique(apply(ped.data.prep$cov.data, 2, class)), "integer")
	expect_equal(levels(ped.data.prep$gen.data[[1]]), c("1", "2", "3", NA))
})

test_that("Preprocessing preprocessed data", {
	haplin.data.prep <- genDataLoad(my.haplin.cov.preproc.out)
	expect_error(genDataPreprocess(data.in = haplin.data.prep, overwrite = if.overwrite))
})

test_that("Preprocessing haplin data, too few markers in map.file", {
	haplin.data.read <- genDataLoad(my.haplin.cov.out)
	my.map <- "tmp.map"
	write.table(
	  matrix(c(1, "rs21", 3), ncol = 3),
	  file = my.map,
	  quote = FALSE,
	  sep = " ",
	  row.names = FALSE,
	  col.names = FALSE
	)
	
	expect_warning(
	  genDataPreprocess(
	    data.in = haplin.data.read,
	    map.file = my.map,
	    overwrite = if.overwrite
	)
	)
})
