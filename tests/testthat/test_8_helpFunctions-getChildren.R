context("Testing helper functions: getChildren, getMothers, getFathers")

source("common_vars.R")

data.haplin.no.cov <- genDataLoad(my.haplin.no.cov.out)
data.haplin.cov <- genDataLoad(my.haplin.cov.out)
data.ped <- genDataLoad(my.ped.out.data)
data.haplin.cov.preproc <- genDataLoad(my.haplin.cov.preproc.out)

# ---- get only children ----
test_that("Getting only children, ped format", {
	children.test <- suppressWarnings(
	  getChildren(
  	  data.ped,
  	  file.out = "my_data_onlyChildren",
  	  dir.out = ".", overwrite = TRUE
  	 )
	)
	# children should have at least one parent in the data
	any.parent <- apply(children.test$cov.data, 1, function(row){
		!(row[ "id.f" ] == "0" & row[ "id.m" ] == "0")
	})
	expect_true(all(any.parent))
})

# currently, only PED format is working with these functions
test_that("Getting only children, haplin format", {
	expect_error(
	  getChildren(
	    data.haplin.cov, file.out = "my_data_onlyChildren", dir.out = ".", overwrite = TRUE
	   )
	 )
})

# ---- get only mothers ----
test_that("Getting only mothers, ped format", {
	mothers.test <- suppressWarnings(
	  getMothers(
	    data.ped, file.out = "my_data_onlyMothers", dir.out = ".", overwrite = TRUE
	  )
	 )
	# parents should not have parents in the data
	no.parents <- apply(mothers.test$cov.data, 1, function(row){
		(row[ "id.f" ] == "0" & row[ "id.m" ] == "0")
	})
	expect_true(all(no.parents))
})

# ---- get only fathers ----
test_that("Getting only fathers, ped format", {
	fathers.test <- suppressWarnings(
	  getFathers(
	    data.ped, file.out = "my_data_onlyFathers", dir.out = ".", overwrite = TRUE
	  )
	 )
	# parents should not have parents in the data
	no.parents <- apply(fathers.test$cov.data, 1, function(row){
		(row[ "id.f" ] == "0" & row[ "id.m" ] == "0")
	})
	expect_true(all(no.parents))
})
