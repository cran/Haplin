context("Testing file reading")

source("common_vars.R")

test_that("Reading in ped data", {
	gen.data.read <- suppressWarnings(
	  genDataRead(file.in = my.ped.data, format = "ped", file.out = my.ped.out.data, overwrite = if.overwrite)
	)

	expect_equal(names(gen.data.read), Haplin:::.haplinEnv$.haplin.data.names)
	expect_true(is.matrix(gen.data.read$cov.data))
	expect_true(is.list(gen.data.read$gen.data))
	expect_true(ff::is.ff(gen.data.read$gen.data[[ 1 ]]))
	expect_s3_class(gen.data.read, "haplin.data")
	
	expect_equal(ncol(gen.data.read$cov.data), 6)
	expect_equal(levels(gen.data.read$gen.data[[1]]), c("G", NA, "C", "T", "A"))
	
	expect_named(gen.data.read$aux, aux.list.names)
	expect_equal(gen.data.read$aux$info$filespecs$format, "ped")
})

test_that("Reading in haplin data without covariates and no n.vars given", {
	expect_error(genDataRead(file.in = my.haplin.data.no.cov, format = "haplin", file.out = my.haplin.no.cov.out, overwrite = if.overwrite))
})

test_that("Reading in haplin data without covariates", {
	gen.data.read <- suppressWarnings(
	  genDataRead(
	    file.in = my.haplin.data.no.cov, format = "haplin", n.vars = 0,
	    file.out = my.haplin.no.cov.out, overwrite = if.overwrite
	  )
	)
	
	expect_equal(names(gen.data.read), Haplin:::.haplinEnv$.haplin.data.names)
	expect_true(is.list(gen.data.read$gen.data))
	expect_true(ff::is.ff(gen.data.read$gen.data[[ 1 ]]))
	expect_s3_class(gen.data.read, "haplin.data")

	expect_null(gen.data.read$cov.data)
	expect_equal(levels(gen.data.read$gen.data[[1]]), c("4", "2", "1", "3", NA, "C", "T"))

	expect_named(gen.data.read$aux, aux.list.names)
	expect_equal(gen.data.read$aux$info$filespecs$format, "haplin")
})

test_that("Reading in haplin data with covariates, and wrong sep given", {
	expect_error(genDataRead(file.in = my.haplin.data.cov, format = "haplin", n.vars = 2, file.out = my.haplin.cov.out, overwrite = if.overwrite))
})
	
test_that("Reading in haplin data with covariates", {
	gen.data.read <- suppressWarnings(
	  genDataRead(
	    file.in = my.haplin.data.cov, format = "haplin", n.vars = 2,
	    allele.sep = "", col.sep = "\t", file.out = my.haplin.cov.out,
	    overwrite = if.overwrite
	  )
	)

	expect_equal(names(gen.data.read), Haplin:::.haplinEnv$.haplin.data.names)
	expect_true(is.matrix(gen.data.read$cov.data))
	expect_true(is.list(gen.data.read$gen.data))
	expect_true(ff::is.ff(gen.data.read$gen.data[[ 1 ]]))
	expect_s3_class(gen.data.read, "haplin.data")

	expect_equal(ncol(gen.data.read$cov.data), 2)
	expect_equal(levels(gen.data.read$gen.data[[1]]), c("A", "G", NA, "B"))

	expect_named(gen.data.read$aux, aux.list.names)
	expect_equal(gen.data.read$aux$info$filespecs$format, "haplin")
})

test_that("Reading in PED data with extra covariate file", {
	gen.data.read <- suppressWarnings(
	  genDataRead(
	    file.in = my.ped.data, format = "ped", file.out = my.ped.out.data,
	    cov.file.in = add.cov.file.name, overwrite = if.overwrite
	  )
	)

	expect_equal(names(gen.data.read), Haplin:::.haplinEnv$.haplin.data.names)
	expect_true(is.data.frame(gen.data.read$cov.data))
	expect_true(is.list(gen.data.read$gen.data))
	expect_true(ff::is.ff(gen.data.read$gen.data[[ 1 ]]))
	expect_s3_class(gen.data.read, "haplin.data")
	
	expect_equal(ncol(gen.data.read$cov.data), 10)
	expect_equal(levels(gen.data.read$gen.data[[1]]), c("G", NA, "C", "T", "A"))
	
	expect_named(gen.data.read$aux, aux.list.names)
	expect_equal(gen.data.read$aux$info$filespecs$format, "ped")
})

test_that("Reading in PED data with extra covariate file, wrong cov.header", {
	expect_error(
	  gen.data.read <- genDataRead(
	    file.in = my.ped.data, format = "ped", file.out = my.ped.out.data,
	    cov.file.in = add.cov.no.head.file.name, cov.header = c("smoke", "vitamin"),
	    overwrite = if.overwrite
	  )
	)
})
