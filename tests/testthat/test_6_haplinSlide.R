context("Testing basic haplinSlide runs")

source("common_vars.R")
set.seed(146)

test_that("Running haplinSlide with simple options, haplin format", {
	skip_on_cran()

	haplin.data.in <- suppressWarnings(genDataLoad(my.haplin.cov.preproc.out))
	cur.args <- args.simple
	
	haplin.res.list <- lapply(cur.args, function(args){
		suppressWarnings(
  		haplinSlide(haplin.data.in, design = args$design, use.missing = args$use.missing,
	  		maternal = args$maternal, reference = args$reference, response = args$response,
		  	verbose = FALSE, printout = FALSE)
		)
	})
	haplin.res.haptables <- lapply(haplin.res.list, haptable)
	
	expect_known_value(
	  haplin.res.haptables,
	  "haplinSlide_test_def-opt_haplin-f.rds",
	  tolerance = tolerance_tests
	)
})

test_that("Running haplinSlide with simple options, ped format", {
	skip_on_cran()

	ped.data.in <- genDataLoad(my.ped.preproc.out)
	cur.args <- args.simple
	
	haplin.res.list <- lapply(cur.args, function(args){
		suppressWarnings(
  		haplinSlide(ped.data.in, design = args$design, use.missing = args$use.missing,
	  		maternal = args$maternal, reference = args$reference, response = args$response,
		  	verbose = FALSE, printout = FALSE)
		)
	})
	haplin.res.haptables <- lapply(haplin.res.list, haptable)
	
	expect_known_value(
	  haplin.res.haptables,
	  "haplinSlide_test_def-opt_ped-f.rds",
	  tolerance = tolerance_tests
	)
})

test_that("Running haplinSlide with xchrom, ped format", {
	skip_on_cran()

	ped.data.in <- genDataLoad(my.ped.preproc.out)
	cur.args <- args.xchrom
	
	haplin.res.list <- lapply(cur.args, function(args){
		suppressWarnings(
		  haplinSlide(ped.data.in, design = args$design, use.missing = args$use.missing,
			  maternal = args$maternal, reference = args$reference, response = args$response,
			  xchrom = args$xchrom, comb.sex = args$comb.sex, verbose = FALSE, printout = FALSE)
		)
	})
	haplin.res.haptables <- lapply(haplin.res.list, haptable)
	
	expect_known_value(
	  haplin.res.haptables,
	  "haplinSlide_test_xchrom_ped-f.rds",
	  tolerance = tolerance_tests
	)
})

test_that("Running haplinSlide with covariates and 'cc' design, ped format", {
	skip_on_cran()

	ped.data.read <- genDataLoad(my.ped.out.data)
	# I need to take only some markers, otherwise it would run for too long!
	ped.data.part <- genDataGetPart(ped.data.read, design = "triad", markers = 110:111,
		overwrite = if.overwrite)
	expect_warning(
  	ped.data.in <- genDataPreprocess(
  	  data.in = ped.data.part,
  	  file.out = my.ped.preproc.cc.out,
  		design = "cc",
  		overwrite = if.overwrite
  	)
	)

	cur.args <- args.cc.cc
	
	haplin.res.list <- lapply(cur.args, function(args){
		suppressWarnings(
  		haplinSlide(ped.data.in, design = args$design, ccvar = args$ccvar,
	  		use.missing = args$use.missing, maternal = args$maternal, reference = args$reference,
		  	response = args$response, verbose = FALSE, printout = FALSE)
		)
	})
	haplin.res.haptables <- lapply(haplin.res.list, haptable)
	
	expect_known_value(
	  haplin.res.haptables,
	  "haplinSlide_test_cc-cc_ped-f.rds",
	  tolerance = tolerance_tests
	)
})
