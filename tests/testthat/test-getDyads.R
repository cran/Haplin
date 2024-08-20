test_that("getDyads works", {
  exmpls.dir <- system.file("extdata", package = "Haplin")
  ped.file <- "qc_all_ok"
  map.file <- file.path(exmpls.dir, "qc_all.map")
  my.ped.data <- genDataLoad(ped.file, dir.in = exmpls.dir)

  expect_warning(
    my.ped.data.dyads <- getDyads(my.ped.data, overwrite = TRUE)
  )
  expect_equal(nindiv(my.ped.data.dyads), 176)
  expect_equal(nfam(my.ped.data.dyads), 88)
})
