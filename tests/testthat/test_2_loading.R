context("Testing loading of .ffData")

source("common_vars.R")

test_that("Loading ped data", {
	gen.data.read <- suppressWarnings(
	  genDataRead(
	    file.in = my.ped.data, format = "ped", file.out = my.ped.out.data,
	    overwrite = if.overwrite
	  )
	)
	
	expect_equal(genDataLoad(my.ped.out.data), gen.data.read)
})

test_that("Loading haplin data, no covariates", {
	gen.data.read <- suppressWarnings(
	  genDataRead(
	    file.in = my.haplin.data.no.cov, format = "haplin", n.vars = 0,
	    file.out = my.haplin.no.cov.out, overwrite = if.overwrite
	  )
	)

	expect_equal(genDataLoad(my.haplin.no.cov.out), gen.data.read)	
})

test_that("Loading haplin data, with covariates", {
	gen.data.read <- suppressWarnings(
	  genDataRead(
	    file.in = my.haplin.data.cov, format = "haplin", n.vars = 2,
	    allele.sep = "", col.sep = "\t", file.out = my.haplin.cov.out,
	    overwrite = if.overwrite
	 )
	)

	expect_equal(genDataLoad(my.haplin.cov.out), gen.data.read)	
})
