
## ---- tissot_raster ----

test_that("tissot_raster returns correct structure", {
  tr <- tissot_raster("+proj=robin", nx = 20)
  expect_s3_class(tr, "tissot_raster")
  expect_length(tr$x, 20L)
  expect_true(length(tr$y) > 0L)
  expect_true(all(c("scale.area", "angle_deformation") %in% names(tr$z)))
  expect_equal(dim(tr$z[[1L]]), c(length(tr$y), length(tr$x)))
})

test_that("tissot_raster with explicit extent", {
  tr <- tissot_raster("+proj=utm +zone=55 +south",
                      extent = c(100000, 900000, 5000000, 6500000),
                      nx = 15)
  expect_equal(length(tr$x), 15L)
  expect_true(any(!is.na(tr$z$scale.area)))
})

test_that("tissot_raster custom metrics", {
  tr <- tissot_raster("+proj=moll", nx = 15,
                      metrics = c("scale.area", "scale.h"))
  expect_equal(names(tr$z), c("scale.area", "scale.h"))
})

test_that("print.tissot_raster works", {
  tr <- tissot_raster("+proj=robin", nx = 10)
  expect_output(print(tr), "Tissot distortion raster")
})

test_that("auto_extent fails gracefully with bad CRS", {
  expect_error(tissot_raster("not_a_crs", nx = 5))
})
