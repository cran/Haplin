context("Testing basic haplinStrat runs")

source("common_vars.R")
set.seed(146)

ped.data.read <- genDataLoad(my.ped.out.data)
strata.col <- sample(x = 1:4, size = nrow(ped.data.read$cov.data), replace = TRUE)
ped.data.read$cov.data <- cbind(ped.data.read$cov.data, strata.col)
colnames(ped.data.read$cov.data)[ ncol(ped.data.read$cov.data) ] <- "strata"

ped.data.part <- suppressWarnings(
  genDataGetPart(
    ped.data.read, design = "triad", markers = 110:111, overwrite = if.overwrite
	)
)
ped.data.in <-suppressWarnings(
  genDataPreprocess(
    data.in = ped.data.part,
    file.out = my.ped.preproc.env.out,
	  overwrite = if.overwrite
   )
 )

test_that("Running haplinStrat with simple options, haplin format", {
	skip_on_cran()

	haplin.data.in <- genDataLoad(my.haplin.cov.preproc.out)
	cur.args <- lapply(args.simple, function(x) { c(x, strata = 1) })

	haplin.res.list <- lapply(cur.args, function(args){
		suppressWarnings(
		  haplinStrat(
  		  haplin.data.in,
        design = args$design,
        markers = 1,
  			use.missing = args$use.missing,
        maternal = args$maternal,
        reference = args$reference,
  			response = args$response,
        strata = args$strata,
        verbose = FALSE,
        printout = FALSE
  		)
		)
	})
	haplin.res.haptables <- lapply(haplin.res.list, haptable)
	
	expect_known_value(
	  haplin.res.haptables,
	  "haplinStrat_test_def-opt_haplin-f.rds",
	  tolerance = tolerance_tests
	 )
})

test_that("Running haplinStrat with simple options, ped format", {
	skip_on_cran()

	cur.args <- lapply(args.simple, function(x) { c(x, strata = 13) })
	
	haplin.res.list <- lapply(cur.args, function(args){
		suppressWarnings(
		  haplinStrat(
  		  ped.data.in,
        design = args$design,
        markers = 1,
  			use.missing = args$use.missing,
        maternal = args$maternal,
        reference = args$reference,
  			response = args$response,
        strata = args$strata,
        verbose = FALSE,
        printout = FALSE
		 )
		)
	})
	haplin.res.haptables <- lapply(haplin.res.list, haptable)
	
	expect_known_value(
	  haplin.res.haptables,
	  "haplinStrat_test_def-opt_ped-f.rds",
	  tolerance = tolerance_tests
	  )
})

test_that("Running haplinStrat with xchrom, ped format", {
	skip_on_cran()

	cur.args <- lapply(args.xchrom, function(x) { c(x, strata = 13) })
	
	haplin.res.list <- lapply(cur.args, function(args){
		suppressWarnings(
		  haplinStrat(
  		  ped.data.in,
        design = args$design,
        markers = 1,
  			use.missing = args$use.missing,
        maternal = args$maternal,
        reference = args$reference,
  			response = args$response,
        xchrom = args$xchrom,
        comb.sex = args$comb.sex,
  			strata = args$strata,
        verbose = FALSE,
        printout = FALSE
  		)
		)
	})
	haplin.res.haptables <- lapply(haplin.res.list, haptable)
	
	expect_known_value(
	  haplin.res.haptables,
	  "haplinStrat_test_xchrom_ped-f.rds",
	  tolerance = tolerance_tests
	  )
})

test_that("Running haplinStrat with covariates and 'cc' design, ped format", {
	skip_on_cran()

	ped.data.in <- suppressWarnings(
	  genDataPreprocess(
	    data.in = ped.data.part, file.out = my.ped.preproc.cc.env.out,
		  design = "cc", overwrite = if.overwrite
		)
	)

	cur.args <- lapply(args.cc.cc, function(x) { c(x, strata = 5) })
	
	haplin.res.list <- lapply(cur.args, function(args){
		suppressWarnings(
		  haplinStrat(
  		  ped.data.in,
        design = args$design,
        markers = 1,
        ccvar = args$ccvar,
  			use.missing = args$use.missing,
        maternal = args$maternal,
        reference = args$reference,
  			response = args$response,
        strata = args$strata,
        verbose = FALSE,
        printout = FALSE
  		)
		)
	})
	haplin.res.haptables <- lapply(haplin.res.list, haptable)
	
	expect_known_value(
	  haplin.res.haptables,
	  "haplinStrat_test_cc-cc_ped-f.rds",
	  tolerance = tolerance_tests
	  )
})
