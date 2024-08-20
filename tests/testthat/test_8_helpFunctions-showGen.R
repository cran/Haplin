context("Testing helper functions: showGen")

source("common_vars.R")
library(ff)

data.haplin.no.cov <- genDataLoad(my.haplin.no.cov.out)
data.haplin.cov <- genDataLoad(my.haplin.cov.out)
data.ped <- genDataLoad(my.ped.out.data)
# # take only 100 individuals, to operate on smaller data
# data.ped <- genDataGetPart(data.ped, design = "triad",
# 	rows = runif(100)*nrow(data.ped$cov.data))
data.haplin.cov.preproc <- genDataLoad(my.haplin.cov.preproc.out)

test_that("Extracting genotypes from preprocessed data", {
	expect_error(showGen(data.haplin.cov.preproc))
})

test_that("Extracting genotypes from haplin-formatted data without covariates", {
	expect_error(showGen(data.haplin.no.cov, sex = 2))
})

# ---- all data ----
test_that("Showing all genotypes from haplin-formatted data", {
	extracted.gen.data <- showGen(data.haplin.cov, markers = "all")
	expected <- as.ff(data.haplin.cov$gen.data[[1]][ 1:5, ], vmode = "byte")

	expect_equal(extracted.gen.data, expected)
})

test_that("Showing all genotypes from ped-formatted data", {
	extracted.gen.data <- showGen(data.ped, markers = "all")
	expected <- as.ff(data.ped$gen.data[[1]][ 1:5, ], vmode = "byte")

	expect_equal(extracted.gen.data, expected)
})

# ---- only males ----
test_that("Extracting genotypes for males, haplin-formatted data", {
	# Haplin format does not have sex of individuals by default
	expect_error(showGen(data.haplin.cov, sex = 1))
})

test_that("Extracting genotypes for males, ped-formatted data", {
	extracted.gen.data <- showGen(data.ped, sex = 1)
	expected <- as.ff(data.ped$gen.data[[1]][ data.ped$cov.data[ ,"sex" ] == "1",1:10 ],
		vmode = "byte")

	expect_equal(extracted.gen.data, expected)
})

test_that("Extracting genotypes for males, ped-formatted data, wrong number given", {
	expect_error(showGen(data.ped, sex = 0))
})

# ---- default: first 5 entries ----
test_that("Showing first 5 entries, haplin format", {
	extracted.gen.data <- showGen(data.haplin.cov, markers = 1:2)
	expected <- as.ff(data.haplin.cov$gen.data[[1]][ 1:5, ], vmode = "byte")
	
	expect_equal(extracted.gen.data, expected)
})

test_that("Showing first 5 entries, ped format", {
	chosen.markers <- c(23:25,101)
	extracted.gen.data <- showGen(data.ped, markers = chosen.markers)
	chosen.columns <- c(sapply(chosen.markers, function(x){
		((x - 1)*2 + 1):(x*2)
	}))
	expected <- as.ff(data.ped$gen.data[[1]][ 1:5,chosen.columns ], vmode = "byte")
	
	expect_equal(extracted.gen.data, expected)
})

