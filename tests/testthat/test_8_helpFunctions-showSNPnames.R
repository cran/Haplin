context("Testing helper functions: showSNPnames")

source("common_vars.R")

data.haplin.no.cov <- genDataLoad(my.haplin.no.cov.out)
data.haplin.cov <- genDataLoad(my.haplin.cov.out)
data.ped <- genDataLoad(my.ped.out.data)
data.haplin.cov.preproc <- genDataLoad(my.haplin.cov.preproc.out)

# ---- all data ----
test_that("Showing all marker names from haplin-formatted data", {
	all.snp.names <- showSNPnames(data.haplin.cov, n = "all")

	expect_equal(all.snp.names, data.haplin.cov$aux$marker.names)
	expect_type(colnames(data.haplin.cov$gen.data[[1]]), "character")
})

test_that("Showing all marker names from ped-formatted data", {
	all.snp.names <- showSNPnames(data.ped, n = "all")

	expect_equal(all.snp.names, data.ped$aux$marker.names)
})

# ---- default: first 5 entries ----
test_that("Showing first entries, haplin format", {
	all.snp.names <- showSNPnames(data.haplin.cov)
	expected <- data.haplin.cov$aux$marker.names[ 1:2 ]
	
	expect_equal(all.snp.names, expected)
})

test_that("Showing first 5 entries, ped format", {
	all.snp.names <- showSNPnames(data.ped)
	expected <- data.ped$aux$marker.names[ 1:5 ]
	
	expect_equal(all.snp.names, expected)
})

# ---- showing first 15 entries, user sets only 'n' ----
test_that("Showing first 2 entries, only n given, haplin format", {
	n <- 2
	all.snp.names <- showSNPnames(data.haplin.cov, n)
	expected <- data.haplin.cov$aux$marker.names[ 1:n ]
	
	expect_equal(all.snp.names, expected)
})

test_that("Showing first 15 entries, only n given, ped format", {
	n <- 15
	all.snp.names <- showSNPnames(data.ped, n)
	expected <- data.ped$aux$marker.names[ 1:n ]
	
	expect_equal(all.snp.names, expected)
})

test_that("Trying to show entries outside the data range, only n given, haplin format", {
	n <- length(data.haplin.cov$aux$marker.names) + 2
	expect_error(showSNPnames(data.haplin.cov, n = n))
})

test_that("Trying to show entries outside the data range, only n given, ped format", {
	n <- length(data.ped$aux$marker.names) + 2
	expect_error(showSNPnames(data.ped, n = n))
})

# ---- showing 15 entries, user sets 'from' and 'n' ----
n <- 15
from <- 10
test_that("Showing 15 entries from 10, only n and from given, ped format", {
	all.snp.names <- showSNPnames(data.ped, n, from = from)
	expected <- data.ped$aux$marker.names[ from:(n + from - 1) ]
	
	expect_equal(all.snp.names, expected)
})

test_that("Trying to show entries outside the data range, only n and from given, haplin format", {
	from <- length(data.haplin.cov$aux$marker.names) - 10
	expect_error(showSNPnames(data.haplin.cov, n = n, from = from))
})

test_that("Trying to show entries outside the data range, only n and from given, ped format", {
	from <- length(data.ped$aux$marker.names) - 10
	expect_error(showSNPnames(data.ped, n = n, from = from))
})

# ---- showing 15 entries, user sets 'to' and 'n' ----
n <- 15
to <- 100
test_that("Showing 15 entries from 10, only to and from given, ped format", {
	all.snp.names <- showSNPnames(data.ped, n, to = to)
	expected <- data.ped$aux$marker.names[ (to - n + 1):to ]
	
	expect_equal(all.snp.names, expected)
})

n <- 115
test_that("Trying to show entries outside the data range, only to and n given, haplin format", {
	expect_error(showSNPnames(data.haplin.cov, n, to = to))
})

test_that("Trying to show entries outside the data range, only to and n given, ped format", {
	expect_error(showSNPnames(data.ped, n, to = to))
})

# ---- showing chosen entries, all parameters given ----
# (should use only 'to' and 'from', and show a warning)
n <- 15
to <- 100
from <- 23
test_that("Showing entries from 23 to 100, all parameters given, ped format", {
	expect_warning(
	  all.snp.names <- showSNPnames(data.ped, n, from = from, to = to)
   )
	expected <- data.ped$aux$marker.names[ from:to ]
	
	expect_equal(all.snp.names, expected)
})

# ---- showing chosen entries, 'to' and 'from' given ----
test_that("Showing entries from 23 to 100, 'to' and 'from' given, ped format", {
	all.snp.names <- showSNPnames(data.ped, from = from, to = to)
	expected <- data.ped$aux$marker.names[ from:to ]
	
	expect_equal(all.snp.names, expected)
})

from <- 0
test_that("Trying to show entries outside the data range, 'to' and 'from' given, haplin format", {
	expect_error(showSNPnames(data.haplin.cov, from = from, to = to))
})

test_that("Trying to show entries outside the data range, 'to' and 'from' given, ped format", {
	expect_error(showSNPnames(data.ped, from = from, to = to))
})

from <- 15
test_that("Trying to show entries outside the data range, 'to' and 'from' given, haplin format", {
	to <- length(data.haplin.cov$aux$marker.names) + 2
	expect_error(showSNPnames(data.haplin.cov, from = from, to = to))
})

test_that("Trying to show entries outside the data range, 'to' and 'from' given, ped format", {
	to <- length(data.ped$aux$marker.names) + 2
	expect_error(showSNPnames(data.ped, from = from, to = to))
})

from <- 25
to <- 13
test_that("'From' parameter larger than 'to', 'to' and 'from' given, haplin format", {
	expect_error(showSNPnames(data.haplin.cov, from = from, to = to))
})

test_that("'From' parameter larger than 'to', 'to' and 'from' given, ped format", {
	expect_error(showSNPnames(data.ped, from = from, to = to))
})
