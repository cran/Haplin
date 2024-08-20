context("Testing basic haplin runs")

source("common_vars.R")
set.seed(146)

test_that("Running haplin with simple options, haplin format", {
	haplin.data.in <- genDataLoad(my.haplin.cov.preproc.out)
	cur.args <- args.simple
	
	haplin.res.list <- lapply(cur.args, function(args){
		suppressWarnings(
		  haplin(haplin.data.in, design = args$design, use.missing = args$use.missing,
			  maternal = args$maternal, reference = args$reference, response = args$response,
			  verbose = FALSE, printout = FALSE)
		)
	})
	haplin.res.haptables <- lapply(haplin.res.list, haptable)
	
	expect_known_value(
	  object = haplin.res.haptables,
	  file = "haplin_test_def-opt_haplin-f.rds",
	  tolerance = tolerance_tests
	 )
})

test_that("Running haplin with simple options, ped format", {
	ped.data.in <- genDataLoad(my.ped.preproc.out)
	cur.args <- args.simple
	
	haplin.res.list <- lapply(cur.args, function(args){
		suppressWarnings(
		  haplin(ped.data.in, design = args$design, use.missing = args$use.missing,
			  maternal = args$maternal, reference = args$reference, response = args$response,
			  verbose = FALSE, printout = FALSE)
		)
	})
	haplin.res.haptables <- lapply(haplin.res.list, haptable)
	
	expect_known_value(
	  haplin.res.haptables,
	  "haplin_test_def-opt_ped-f.rds",
	  tolerance = tolerance_tests
	 )
})

test_that("Running haplin with xchrom, ped format", {
	ped.data.in <- genDataLoad(my.ped.preproc.out)
	cur.args <- args.xchrom
	
	haplin.res.list <- lapply(cur.args, function(args){
		suppressWarnings(
		  haplin(ped.data.in, design = args$design, use.missing = args$use.missing,
			  maternal = args$maternal, reference = args$reference, response = args$response,
			  xchrom = args$xchrom, comb.sex = args$comb.sex, verbose = FALSE, printout = FALSE)
		)
	})
	haplin.res.haptables <- lapply(haplin.res.list, haptable)
	
	expect_known_value(
	  haplin.res.haptables,
	  "haplin_test_xchrom_ped-f.rds",
	  tolerance = tolerance_tests
	 )
})

test_that("Running haplin with covariates and 'cc' design, ped format", {
	ped.data.read <- genDataLoad(my.ped.out.data)
	# I need to take only some markers, otherwise it would run for too long!
	ped.data.part <- genDataGetPart(ped.data.read, design = "triad", markers = 110:111,
		overwrite = if.overwrite)
	ped.data.in <- suppressWarnings(
		  genDataPreprocess(
		    data.in = ped.data.part, file.out = my.ped.preproc.cc.out, design = "cc",
		    overwrite = if.overwrite
		 )
		)

	cur.args <- args.cc.cc
	
	haplin.res.list <- lapply(cur.args, function(args){
		suppressWarnings(
  		haplin(ped.data.in, design = args$design, ccvar = args$ccvar,
			  use.missing = args$use.missing, maternal = args$maternal, reference = args$reference,
			  response = args$response, verbose = FALSE, printout = FALSE)
		)
	})
	haplin.res.haptables <- lapply(haplin.res.list, haptable)
	
	expect_known_value(
	  haplin.res.haptables,
	  "haplin_test_cc-cc_ped-f.rds",
	  tolerance = tolerance_tests
	 )
})
