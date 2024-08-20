context("Testing helper functions: getFullTriads and getDyads")

source("common_vars.R")

data.haplin.no.cov <- genDataLoad(my.haplin.no.cov.out)
data.haplin.cov <- genDataLoad(my.haplin.cov.out)
data.ped <- genDataLoad(my.ped.out.data)
data.haplin.cov.preproc <- genDataLoad(my.haplin.cov.preproc.out)

# ---- get all triads ----
test_that("Getting only full triads, ped format", {
	triads.test <- suppressWarnings(
	  getFullTriads(
	    data.ped, file.out = "my_data_onlyTriads", dir.out = ".", overwrite = TRUE
	 )
	)
	
	# the number of individuals in the data should be divisible by 3
	expect_true(nindiv(triads.test) %% 3 == 0)
	# there should not be any row with only NAs in the genetic data
	find.NAs <- lapply(triads.test$gen.data, function(gen.el){
			which.na.level <- which(is.na(levels(gen.el)))
			apply(gen.el[,], 1, function(row){
				# how many are NOT NAs
				sum(as.numeric(row) != which.na.level)
			})
	})
	if(length(find.NAs) == 1){
		sum.not.na <- matrix(find.NAs[[1]], ncol = 1)
	} else {
		sum.not.na <- Reduce(cbind, find.NAs)
		colnames(sum.not.na) <- NULL
	}
	# this tells us which rows have no NAs
	rows.only.nas <- rowSums(sum.not.na) == 0
	expect_true(all(!rows.only.nas))
})

# ---- get only dyads ----
test_that("Getting only dyads, ped format", {
	dyads.test <- suppressWarnings(
	  getDyads(
	    data.ped, file.out = "my_data_onlyDyads", dir.out = ".", overwrite = TRUE
	  )
	)
	
	# the number of individuals in the data should be divisible by 2
	expect_true(nindiv(dyads.test) %% 2 == 0)
	# there should not be any row with only NAs in the genetic data
	find.NAs <- lapply(dyads.test$gen.data, function(gen.el){
			which.na.level <- which(is.na(levels(gen.el)))
			apply(gen.el[,], 1, function(row){
				# how many are NOT NAs
				sum(as.numeric(row) != which.na.level)
			})
	})
	if(length(find.NAs) == 1){
		sum.not.na <- matrix(find.NAs[[1]], ncol = 1)
	} else {
		sum.not.na <- Reduce(cbind, find.NAs)
		colnames(sum.not.na) <- NULL
	}
	# this tells us which rows have no NAs
	rows.only.nas <- rowSums(sum.not.na) == 0
	expect_true(all(!rows.only.nas))
})
