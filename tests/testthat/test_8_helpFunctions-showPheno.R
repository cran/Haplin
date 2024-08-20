context("Testing helper functions: showPheno")

source("common_vars.R")

data.haplin.no.cov <- genDataLoad(my.haplin.no.cov.out)
data.haplin.cov <- genDataLoad(my.haplin.cov.out)
data.ped <- genDataLoad(my.ped.out.data)
data.haplin.cov.preproc <- genDataLoad(my.haplin.cov.preproc.out)

test_that("Extracting phenotype data from preprocessed data", {
	expect_error(showPheno(data.haplin.cov.preproc))
})

test_that("Extracting phenotype data from haplin-formatted data without covariates", {
	expect_error(showPheno(data.haplin.no.cov))
	expect_error(showPheno(data.haplin.no.cov, n = 2))
})

# ---- all data ----
test_that("Showing all phenotype data from haplin-formatted data", {
	extracted.cov.data <- showPheno(data.haplin.cov, n = "all")

	expect_equal(extracted.cov.data, data.haplin.cov$cov.data)
})

test_that("Showing all phenotype data from ped-formatted data", {
	extracted.cov.data <- showPheno(data.ped, n = "all")

	expect_equal(extracted.cov.data, data.ped$cov.data)
})

# ---- only males ----
test_that("Extracting covariate data for males, haplin-formatted data", {
	# Haplin format does not have sex of individuals by default
	expect_error(showPheno(data.haplin.cov, sex = 1))
})

test_that("Extracting covariate data for males, ped-formatted data", {
	extracted.cov.data <- showPheno(data.ped, sex = 1)
	expected <- data.ped$cov.data[ data.ped$cov.data[ ,"sex" ] == 1, ]

	expect_equal(extracted.cov.data, expected)
})

test_that("Extracting covariate data for males, ped-formatted data, wrong number given", {
	expect_error(showPheno(data.ped, sex = 0))
})

# ---- default: first 5 entries ----
test_that("Showing first 5 entries, haplin format", {
	extracted.cov.data <- showPheno(data.haplin.cov)
	expected <- data.haplin.cov$cov.data[ 1:5, ]
	
	expect_equal(extracted.cov.data, expected)
})

test_that("Showing first 5 entries, ped format", {
	extracted.cov.data <- showPheno(data.ped)
	expected <- data.ped$cov.data[ 1:5, ]
	
	expect_equal(extracted.cov.data, expected)
})

# ---- showing first 15 entries, user sets only 'n' ----
n <- 15
test_that("Showing first 15 entries, only n given, haplin format", {
	extracted.cov.data <- showPheno(data.haplin.cov, n)
	expected <- data.haplin.cov$cov.data[ 1:n, ]
	
	expect_equal(extracted.cov.data, expected)
})

test_that("Showing first 15 entries, only n given, ped format", {
	extracted.cov.data <- showPheno(data.ped, n)
	expected <- data.ped$cov.data[ 1:n, ]
	
	expect_equal(extracted.cov.data, expected)
})

test_that("Trying to show entries outside the data range, only n given, haplin format", {
	n <- nrow(data.haplin.cov$cov.data) + 2
	expect_error(showPheno(data.haplin.cov, n))
})

test_that("Trying to show entries outside the data range, only n given, ped format", {
	n <- nrow(data.ped$cov.data) + 2
	expect_error(showPheno(data.ped, n, from = from))
})

# ---- showing 15 entries, user sets 'from' and 'n' ----
n <- 15
from <- 10
test_that("Showing 15 entries from 10, only n and from given, haplin format", {
	extracted.cov.data <- showPheno(data.haplin.cov, n, from = from)
	expected <- data.haplin.cov$cov.data[ from:(n + from - 1), ]
	
	expect_equal(extracted.cov.data, expected)
})

test_that("Showing 15 entries from 10, only n and from given, ped format", {
	extracted.cov.data <- showPheno(data.ped, n, from = from)
	expected <- data.ped$cov.data[ from:(n + from - 1), ]
	
	expect_equal(extracted.cov.data, expected)
})

test_that("Trying to show entries outside the data range, only n and from given, haplin format", {
	from <- nrow(data.haplin.cov$cov.data) - 10
	expect_error(showPheno(data.haplin.cov, n, from = from))
})

test_that("Trying to show entries outside the data range, only n and from given, ped format", {
	from <- nrow(data.ped$cov.data) - 10
	expect_error(showPheno(data.ped, n, from = from))
})

# ---- showing 15 entries, user sets 'to' and 'n' ----
n <- 15
to <- 100
test_that("Showing 15 entries from 10, only to and from given, haplin format", {
	extracted.cov.data <- showPheno(data.haplin.cov, n, to = to)
	expected <- data.haplin.cov$cov.data[ (to - n + 1):to, ]
	
	expect_equal(extracted.cov.data, expected)
})

test_that("Showing 15 entries from 10, only to and from given, ped format", {
	extracted.cov.data <- showPheno(data.ped, n, to = to)
	expected <- data.ped$cov.data[ (to - n + 1):to, ]
	
	expect_equal(extracted.cov.data, expected)
})

n <- 115
test_that("Trying to show entries outside the data range, only to and n given, haplin format", {
	expect_error(showPheno(data.haplin.cov, n, to = to))
})

test_that("Trying to show entries outside the data range, only to and n given, ped format", {
	expect_error(showPheno(data.ped, n, to = to))
})

# ---- showing chosen entries, all parameters given ----
# (should use only 'to' and 'from', and show a warning)
n <- 15
to <- 100
from <- 23
test_that("Showing entries from 23 to 100, all parameters given, haplin format", {
	expect_warning(
	  extracted.cov.data <- showPheno(data.haplin.cov, n, from = from, to = to)
  )
	expected <- data.haplin.cov$cov.data[ from:to, ]
	
	expect_equal(extracted.cov.data, expected)
})

test_that("Showing entries from 23 to 100, all parameters given, ped format", {
	expect_warning(
	 extracted.cov.data <- showPheno(data.ped, n, from = from, to = to)
	)
	expected <- data.ped$cov.data[ from:to, ]
	
	expect_equal(extracted.cov.data, expected)
})

# ---- showing chosen entries, 'to' and 'from' given ----
test_that("Showing entries from 23 to 100, 'to' and 'from' given, haplin format", {
	extracted.cov.data <- showPheno(data.haplin.cov, from = from, to = to)
	expected <- data.haplin.cov$cov.data[ from:to, ]
	
	expect_equal(extracted.cov.data, expected)
})

test_that("Showing entries from 23 to 100, 'to' and 'from' given, ped format", {
	extracted.cov.data <- showPheno(data.ped, from = from, to = to)
	expected <- data.ped$cov.data[ from:to, ]
	
	expect_equal(extracted.cov.data, expected)
})

from <- 0
test_that("Trying to show entries outside the data range, 'to' and 'from' given, haplin format", {
	expect_error(showPheno(data.haplin.cov, from = from, to = to))
})

test_that("Trying to show entries outside the data range, 'to' and 'from' given, ped format", {
	expect_error(showPheno(data.ped, from = from, to = to))
})

from <- 15
test_that("Trying to show entries outside the data range, 'to' and 'from' given, haplin format", {
	to <- nrow(data.haplin.cov$cov.data) + 2
	expect_error(showPheno(data.haplin.cov, from = from, to = to))
})

test_that("Trying to show entries outside the data range, 'to' and 'from' given, ped format", {
	to <- nrow(data.ped$cov.data) + 2
	expect_error(showPheno(data.ped, from = from, to = to))
})

from <- 25
to <- 13
test_that("'From' parameter larger than 'to', 'to' and 'from' given, haplin format", {
	expect_error(showPheno(data.haplin.cov, from = from, to = to))
})

test_that("'From' parameter larger than 'to', 'to' and 'from' given, ped format", {
	expect_error(showPheno(data.ped, from = from, to = to))
})

