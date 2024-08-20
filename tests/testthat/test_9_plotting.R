context("Testing plotting")

source("common_vars.R")
set.seed(146)

ped.data.in <- genDataLoad(my.ped.preproc.out)
cur.args <- args.simple

haplin.res.list <- lapply(cur.args, function(args){
	suppressWarnings(
	  haplinSlide(
	    ped.data.in, design = args$design, use.missing = args$use.missing,
		  maternal = args$maternal, reference = args$reference, response = args$response,
		  verbose = FALSE, printout = FALSE, table.output = TRUE
		)
	)
})

haplin.no.tab.res.list <- lapply(cur.args, function(args){
	suppressWarnings(
	  haplinSlide(
	    ped.data.in, design = args$design, maternal = args$maternal,
	    reference = args$reference, response = args$response,
		  use.missing = args$use.missing, verbose = FALSE, printout = FALSE,
		  table.output = FALSE
		)
	)
})

test_that("Plotting haplinSlide, generated with table.output = TRUE", {
	skip_on_cran()

	# checking that there were no errors
	hapSlide_plots <- lapply(haplin.res.list, plot)
	lapply(hapSlide_plots, function(x){
	  expect_s3_class(x, "ggplot")
	})
	expect_known_value(
	  hapSlide_plots,
	  "hapSlide_plots_table_out_simple.rds",
	  tolerance = tolerance_tests
	  )
})

test_that("Plotting haplinSlide, generated with table.output = FALSE", {
	skip_on_cran()

	# checking that there were no errors
	hapSlide_plots <- lapply(haplin.no.tab.res.list, plot)
	lapply(hapSlide_plots, function(x){
	  expect_s3_class(x, "ggplot")
	})

	expect_known_value(
	  hapSlide_plots,
		"hapSlide_plots_no_table_out_simple.rds",
	  tolerance = tolerance_tests
	)
})

test_that("Plotting haplinSlide, generated with xchrom option, ped format, table.output = TRUE", {
	skip_on_cran()

	cur.args <- args.xchrom
	
	haplin.res.list <- lapply(cur.args, function(args){
		suppressWarnings(
		  haplinSlide(
		    ped.data.in, design = args$design, use.missing = args$use.missing,
			  maternal = args$maternal, reference = args$reference, response = "mult",
			  xchrom = args$xchrom, comb.sex = args$comb.sex, verbose = FALSE,
			  printout = FALSE
			)
		)
	})

	# checking that there were no errors
	hapSlide_plots <- lapply(haplin.res.list, plot)
	lapply(hapSlide_plots, function(x){
	  expect_s3_class(x, "ggplot")
	})
	expect_known_value(
	  hapSlide_plots,
		"hapSlide_plots_table_out_xchrom.rds",
	  tolerance = tolerance_tests
		)
})

test_that("Plotting haplinSlide, generated with covariates and 'cc' design, ped format, table.output = TRUE", {
	skip_on_cran()

	ped.data.read <- genDataLoad(my.ped.out.data)
	# I need to take only some markers, otherwise it would run for too long!
	ped.data.part <- suppressWarnings(
	  genDataGetPart(
	    ped.data.read, design = "triad", markers = 110:111, overwrite = if.overwrite
	   )
	)
	ped.data.in <- suppressWarnings(
	  genDataPreprocess(
	    data.in = ped.data.part, file.out = my.ped.preproc.cc.out,
		  design = "cc", overwrite = if.overwrite
		)
	)

	cur.args <- args.cc.cc
	
	haplin.res.list <- lapply(cur.args, function(args){
		suppressWarnings(
		  haplinSlide(
		    ped.data.in, design = args$design, ccvar = args$ccvar, 
			  use.missing = args$use.missing, maternal = args$maternal,
			  reference = args$reference, response = args$response,
			  verbose = FALSE, printout = FALSE
			 )
		)
	})

	# checking that there were no errors
	hapSlide_plots <- lapply(haplin.res.list, plot)
	lapply(hapSlide_plots, function(x){
	  expect_s3_class(x, "ggplot")
	})
	expect_known_value(
	  hapSlide_plots,
		"hapSlide_plots_table_out_cc.rds",
	  tolerance = tolerance_tests
		)
})
