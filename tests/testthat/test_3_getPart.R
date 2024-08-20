context("Testing fetching parts of data")

source("common_vars.R")
# this is only to produce the .ffData and .RData files
gen.data.read <- suppressWarnings(
  genDataRead(
    file.in = my.ped.data, format = "ped", file.out = my.ped.out.data,
    overwrite = if.overwrite
 )
)
gen.data.read <- suppressWarnings(
  genDataRead(
    file.in = my.haplin.data.no.cov, format = "haplin", n.vars = 0,
    file.out = my.haplin.no.cov.out, overwrite = if.overwrite
 )
)
rm(gen.data.read)

# ---- PED data ----
test_that("Getting markers from PED data", {
	ped.data.read <- suppressWarnings(genDataLoad(my.ped.out.data))
	ped.data.part <- suppressWarnings(
	  genDataGetPart(ped.data.read, design = "triad", markers = 31:40,
	                  overwrite = if.overwrite)
	)
	
	expect_s3_class(ped.data.part, "haplin.data")
	expect_equal(ncol(ped.data.part$gen.data[[1]]), 20)
	expect_equal(nrow(ped.data.part$gen.data[[1]]), nrow(ped.data.read$gen.data[[1]]))
})

test_that("Getting markers from PED data, outside range", {
	ped.data.read <- suppressWarnings(genDataLoad(my.ped.out.data))
	expect_error(genDataGetPart(ped.data.read, design = "triad", markers = 425:430, overwrite = if.overwrite))
})

test_that("Getting individuals from PED data", {
	ped.data.read <- suppressWarnings(genDataLoad(my.ped.out.data))
	ped.data.part <- suppressWarnings(
	  genDataGetPart(
	    ped.data.read, design = "triad", indiv.ids = "1", overwrite = if.overwrite
	  )
	)
	
	expect_s3_class(ped.data.part, "haplin.data")
	expect_equal(ncol(ped.data.part$gen.data[[1]]), ncol(ped.data.read$gen.data[[1]]))
	expect_equal(nrow(ped.data.part$gen.data[[1]]), 550)
})

test_that("Getting individuals from PED data, wrong ID", {
	ped.data.read <- suppressWarnings(genDataLoad(my.ped.out.data))
	expect_error(genDataGetPart(ped.data.read, design = "triad", indiv.ids = "100", overwrite = if.overwrite))
})

test_that("Getting rows from PED data", {
	ped.data.read <- suppressWarnings(genDataLoad(my.ped.out.data))
	ped.data.part <- suppressWarnings(
	  genDataGetPart(
	    ped.data.read, design = "triad", rows = 2:10, overwrite = if.overwrite
	  )
	)
	
	expect_s3_class(ped.data.part, "haplin.data")
	expect_equal(ncol(ped.data.part$gen.data[[1]]), ncol(ped.data.read$gen.data[[1]]))
	expect_equal(nrow(ped.data.part$gen.data[[1]]), 9)
})

test_that("Getting cases from PED data", {
	ped.data.read <- suppressWarnings(genDataLoad(my.ped.out.data))
	ped.data.part <- suppressWarnings(
	  genDataGetPart(
	    ped.data.read, design = "triad", cc = "1", overwrite = if.overwrite
	  )
	)
	
	expect_s3_class(ped.data.part, "haplin.data")
	expect_equal(ncol(ped.data.part$gen.data[[1]]), ncol(ped.data.read$gen.data[[1]]))
	expect_equal(nrow(ped.data.part$gen.data[[1]]), 815)
})

test_that("Getting cases from PED data, wrong cc given", {
	ped.data.read <- suppressWarnings(genDataLoad(my.ped.out.data))
	expect_error(genDataGetPart(ped.data.read, design = "triad", cc = "2", overwrite = if.overwrite))
})

test_that("Getting males from PED data", {
	ped.data.read <- suppressWarnings(genDataLoad(my.ped.out.data))
	ped.data.part <- suppressWarnings(
	  genDataGetPart(
	    ped.data.read, design = "triad", sex = "1", overwrite = if.overwrite
	)
	)
	
	expect_s3_class(ped.data.part, "haplin.data")
	expect_equal(ncol(ped.data.part$gen.data[[1]]), ncol(ped.data.read$gen.data[[1]]))
	expect_equal(nrow(ped.data.part$gen.data[[1]]), 794)
})

test_that("Getting males from PED data, wrong sex given", {
	ped.data.read <- suppressWarnings(genDataLoad(my.ped.out.data))
	expect_error(genDataGetPart(ped.data.read, design = "triad", sex = "3", overwrite = if.overwrite))
})

test_that("Getting males and cases from PED data", {
	ped.data.read <- suppressWarnings(genDataLoad(my.ped.out.data))
	ped.data.part <- suppressWarnings(
	  genDataGetPart(
	    ped.data.read, design = "triad", sex = "1", cc = "1",
	    overwrite = if.overwrite
	  )
	)
	
	expect_s3_class(ped.data.part, "haplin.data")
	expect_equal(ncol(ped.data.part$gen.data[[1]]), ncol(ped.data.read$gen.data[[1]]))
	expect_equal(nrow(ped.data.part$gen.data[[1]]), 384)
})

test_that("Getting specific markers for males and cases from PED data", {
	ped.data.read <- suppressWarnings(genDataLoad(my.ped.out.data))
	ped.data.part <- suppressWarnings(
	  genDataGetPart(
	    ped.data.read, design = "triad", markers = 54:56, sex = "1", cc = "1",
	    overwrite = if.overwrite
	  )
	)
	
	expect_s3_class(ped.data.part, "haplin.data")
	expect_equal(ncol(ped.data.part$gen.data[[1]]), 6)
	expect_equal(nrow(ped.data.part$gen.data[[1]]), 384)
})

# ---- haplin data ----
gen.data.read <- suppressWarnings(
  genDataRead(
    file.in = my.haplin.data.no.cov, format = "haplin", n.vars = 0,
    file.out = my.haplin.no.cov.out, overwrite = if.overwrite
 )
)
gen.data.read <- suppressWarnings(
  genDataRead(
    file.in = my.haplin.data.cov, format = "haplin", n.vars = 2, allele.sep = "",
    col.sep = "\t", file.out = my.haplin.cov.out, overwrite = if.overwrite
 )
)
rm(gen.data.read)

test_that("Getting markers from haplin data", {
	haplin.data.read <- suppressWarnings(genDataLoad(my.haplin.cov.out))
	haplin.data.part <- suppressWarnings(
	  genDataGetPart(
	    haplin.data.read, design = "triad", markers = 2, overwrite = if.overwrite
	  )
	)
	
	expect_s3_class(haplin.data.part, "haplin.data")
	expect_equal(ncol(haplin.data.part$gen.data[[1]]), 6)
	expect_equal(nrow(haplin.data.part$gen.data[[1]]), nrow(haplin.data.read$gen.data[[1]]))
})

test_that("Getting markers from haplin data, outside range", {
	haplin.data.read <- suppressWarnings(genDataLoad(my.haplin.cov.out))
	expect_error(genDataGetPart(haplin.data.read, design = "triad", markers = 3:4, overwrite = if.overwrite))
})

test_that("Getting individuals from haplin data", {
	haplin.data.read <- suppressWarnings(genDataLoad(my.haplin.cov.out))
	expect_error(genDataGetPart(haplin.data.read, design = "triad", indiv.ids = "1", overwrite = if.overwrite))
})

test_that("Getting rows from haplin data", {
	haplin.data.read <- suppressWarnings(genDataLoad(my.haplin.cov.out))
	haplin.data.part <- suppressWarnings(
	  genDataGetPart(
	    haplin.data.read, design = "triad", rows = 2:10, overwrite = if.overwrite
	  )
	)
	
	expect_s3_class(haplin.data.part, "haplin.data")
	expect_equal(ncol(haplin.data.part$gen.data[[1]]), ncol(haplin.data.read$gen.data[[1]]))
	expect_equal(nrow(haplin.data.part$gen.data[[1]]), 9)
})

test_that("Getting cases from haplin data", {
	haplin.data.read <- suppressWarnings(genDataLoad(my.haplin.cov.out))
	expect_error(genDataGetPart(haplin.data.read, design = "triad", cc = "1", overwrite = if.overwrite))
})

test_that("Getting based on cov.1 from haplin data", {
	haplin.data.read <- suppressWarnings(genDataLoad(my.haplin.cov.out))
	haplin.data.part <- suppressWarnings(
	  genDataGetPart(
	    haplin.data.read, design = "triad", cov.1 = c(1,3), overwrite = if.overwrite
	  )
	)
	
	expect_s3_class(haplin.data.part, "haplin.data")
	expect_equal(ncol(haplin.data.part$gen.data[[1]]), ncol(haplin.data.read$gen.data[[1]]))
	expect_equal(nrow(haplin.data.part$gen.data[[1]]), 286)
})

test_that("Getting based on cov.2 from haplin data", {
	haplin.data.read <- suppressWarnings(genDataLoad(my.haplin.cov.out))
	haplin.data.part <- suppressWarnings(
	  genDataGetPart(
	    haplin.data.read, design = "triad", cov.2 = 1, overwrite = if.overwrite
	  )
	)
	
	expect_s3_class(haplin.data.part, "haplin.data")
	expect_equal(ncol(haplin.data.part$gen.data[[1]]), ncol(haplin.data.read$gen.data[[1]]))
	expect_equal(nrow(haplin.data.part$gen.data[[1]]), 377)
})

test_that("Getting based on cov.2 from haplin data, wrong value", {
	haplin.data.read <- suppressWarnings(genDataLoad(my.haplin.cov.out))
	expect_error(genDataGetPart(haplin.data.read, design = "triad", cov.2 = 3, overwrite = if.overwrite))
})

test_that("Getting based on cov.1 and cov.2 from haplin data", {
	haplin.data.read <- suppressWarnings(genDataLoad(my.haplin.cov.out))
	haplin.data.part <- suppressWarnings(
	  genDataGetPart(
	    haplin.data.read, design = "triad", cov.1 = c(1,3), cov.2 = 0,
	    overwrite = if.overwrite
	  )
	)
	
	expect_s3_class(haplin.data.part, "haplin.data")
	expect_equal(ncol(haplin.data.part$gen.data[[1]]), ncol(haplin.data.read$gen.data[[1]]))
	expect_equal(nrow(haplin.data.part$gen.data[[1]]), 179)
})

test_that("Getting specific markers for selected rows from haplin data", {
	haplin.data.read <- suppressWarnings(genDataLoad(my.haplin.cov.out))
	haplin.data.part <- suppressWarnings(
	  genDataGetPart(
	    haplin.data.read, design = "triad", markers = 1, cov.1 = c(1,3),
	    cov.2 = 0, overwrite = if.overwrite
	  )
	)
	
	expect_s3_class(haplin.data.part, "haplin.data")
	expect_equal(ncol(haplin.data.part$gen.data[[1]]), 6)
	expect_equal(nrow(haplin.data.part$gen.data[[1]]), 179)
})

test_that("Getting based on cov.1 from haplin data, no covariates", {
	haplin.data.read <- suppressWarnings(genDataLoad(my.haplin.no.cov.out))
	expect_error(genDataGetPart(haplin.data.read, design = "triad", cov.1 = 3, overwrite = if.overwrite))
})
