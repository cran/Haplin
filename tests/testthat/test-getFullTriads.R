test_that("getFullTriads works", {
  exmpls.dir <- system.file("extdata", package = "Haplin")
  ped.file <- "qc_all_ok"
  map.file <- file.path(exmpls.dir, "qc_all.map")
  my.ped.data <- suppressWarnings(genDataLoad(ped.file, dir.in = exmpls.dir))

  expect_warning(
    my.ped.data.triads <- getFullTriads(
      my.ped.data,
      overwrite = TRUE
    )
  )
  expect_equal(nindiv(my.ped.data.triads), 1362)
  expect_equal(nfam(my.ped.data.triads), 454)
})
