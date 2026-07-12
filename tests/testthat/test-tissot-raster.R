
## ---- tissot_raster ----

test_that("tissot_raster returns correct structure", {
  tr <- tissot_raster("+proj=robin", nx = 20)
  expect_s3_class(tr, "tissot_raster")
  expect_named(tr, c("x", "y", "z", "target", "extent"))
  expect_length(tr$x, 20L)
  expect_true(length(tr$y) > 0L)
  expect_true(all(c("scale_area", "angle_deformation", "scale_h", "scale_k")
                  %in% names(tr$z)))
  expect_equal(dim(tr$z[[1L]]), c(length(tr$y), length(tr$x)))
})

test_that("tissot_raster default metrics produce non-NA values", {
  tr <- tissot_raster("+proj=robin", nx = 20)
  for (m in names(tr$z)) {
    expect_true(any(is.finite(tr$z[[m]])),
                info = paste("all NA for metric:", m))
  }
})

test_that("tissot_raster with explicit extent", {
  tr <- tissot_raster("+proj=utm +zone=55 +south",
                      extent = c(100000, 900000, 5000000, 6500000),
                      nx = 15)
  expect_equal(length(tr$x), 15L)
  expect_true(any(!is.na(tr$z$scale_area)))
})

test_that("tissot_raster custom metrics", {
  tr <- tissot_raster("+proj=moll", nx = 15,
                      metrics = c("scale_area", "scale_h"))
  expect_equal(names(tr$z), c("scale_area", "scale_h"))
  expect_true(any(is.finite(tr$z[["scale_area"]])))
})

test_that("tissot_raster convergence metric works", {
  tr <- tissot_raster("+proj=moll", nx = 15, metrics = "convergence")
  expect_equal(names(tr$z), "convergence")
  expect_true(any(is.finite(tr$z[["convergence"]])))
})

test_that("print.tissot_raster works", {
  tr <- tissot_raster("+proj=robin", nx = 10)
  expect_output(print(tr), "Tissot distortion raster")
})

test_that("auto_extent fails gracefully with bad CRS", {
  expect_error(tissot_raster("not_a_crs", nx = 5))
})
