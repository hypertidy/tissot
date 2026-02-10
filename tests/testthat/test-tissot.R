## ---- tissot() core engine ----

test_that("conformal projection has near-zero angular deformation", {
  ## Mercator is conformal: scale.a ≈ scale.b, angle_deformation ≈ 0
  r <- tissot(cbind(0, 45), "+proj=merc")
  expect_lt(r$angle_deformation, 0.01)
  expect_equal(r$scale.a, r$scale.b, tolerance = 1e-4)
})

test_that("equal-area projection has areal scale ≈ 1", {
  ## Mollweide is equal-area
  r <- tissot(cbind(0, 45), "+proj=moll")
  expect_equal(r$scale.area, 1.0, tolerance = 5e-3)
})

test_that("Plate Carrée at equator: scale.h ≈ 1, scale.k ≈ 1", {
  r <- tissot(cbind(0, 0), "+proj=eqc")
  expect_equal(r$scale.h, 1.0, tolerance = 0.01)
  expect_equal(r$scale.k, 1.0, tolerance = 0.01)
})

test_that("vectorized: multiple points return correct nrow", {
  xy <- expand.grid(seq(-150, 150, by = 30), seq(-60, 60, by = 30))
  r <- tissot(xy, "+proj=robin")
  expect_equal(nrow(r), nrow(xy))
})

test_that("matrix input == expand.grid input", {
  m <- cbind(c(-30, 0, 30), c(-45, 0, 45))
  r1 <- tissot(m, "+proj=robin")
  r2 <- tissot(as.data.frame(m), "+proj=robin")
  expect_equal(r1$scale.area, r2$scale.area, tolerance = 1e-10)
})

test_that("single point via length-2 vector", {
  r <- tissot(c(147, -42.5), "+proj=utm +zone=55 +south")
  expect_equal(nrow(r), 1L)
  expect_true(all(is.finite(unlist(r[, 7:14]))))
})

test_that("equal-area batch: all scale.area ≈ 1", {
  xy <- expand.grid(seq(-120, 120, by = 60), seq(-60, 60, by = 30))
  r <- tissot(xy, "+proj=cea")
  expect_equal(r$scale.area, rep(1, nrow(xy)), tolerance = 5e-3)
})

test_that("polar points don't produce NaN", {
  r <- tissot(cbind(0, 89), "+proj=stere +lat_0=90")
  expect_true(all(is.finite(r$scale.a)))
})

test_that("custom source CRS works", {
  r <- tissot(cbind(0, 45), "+proj=robin", source = "EPSG:4267")
  expect_equal(nrow(r), 1L)
  expect_true(is.finite(r$scale.area))
})

## ---- tissot_tbl class ----

test_that("tissot_tbl has correct class and attributes", {
  r <- tissot(cbind(0, 45), "+proj=robin")
  expect_s3_class(r, "tissot_tbl")
  expect_equal(attr(r, "source"), "EPSG:4326")
  expect_equal(attr(r, "target"), "+proj=robin")
})

test_that("print.tissot_tbl runs without error", {
  r <- tissot(cbind(0, 45), "+proj=robin")
  expect_output(print(r), "Tissot indicatrix")
})

test_that("summary.tissot_tbl runs without error", {
  r <- tissot(expand.grid(seq(-90, 90, by = 45), seq(-45, 45, by = 45)),
              "+proj=robin")
  expect_output(summary(r), "Target CRS")
})

## ---- indicatrix() ----

test_that("indicatrix from tissot_tbl extracts target", {
  r <- tissot(cbind(0, 45), "+proj=robin")
  ii <- indicatrix(r)
  expect_s3_class(ii, "indicatrix_list")
  expect_equal(attr(ii, "target"), "+proj=robin")
  expect_length(ii, 1L)
})

test_that("indicatrix from raw coords requires target", {
  expect_error(indicatrix(cbind(0, 45)), "target is required")
})

test_that("indicatrix from raw coords works", {
  ii <- indicatrix(cbind(0, 45), "+proj=stere +lat_0=90")
  expect_s3_class(ii, "indicatrix_list")
  expect_length(ii, 1L)
})

## ---- indicatrix_list methods ----

test_that("subsetting preserves class and attributes", {
  xy <- expand.grid(seq(-90, 90, by = 45), seq(-45, 45, by = 45))
  ii <- indicatrix(xy, "+proj=robin")
  sub <- ii[1:3]
  expect_s3_class(sub, "indicatrix_list")
  expect_equal(attr(sub, "target"), "+proj=robin")
  expect_length(sub, 3L)
})

test_that("print.indicatrix_list works", {
  ii <- indicatrix(cbind(0, 45), "+proj=robin")
  expect_output(print(ii), "Indicatrix list")
})

## ---- ti_ellipse ----

test_that("ti_ellipse returns closed polygon", {
  ii <- indicatrix(cbind(0, 45), "+proj=robin")
  ell <- ti_ellipse(ii[[1]], scale = 1e5, n = 36)
  expect_equal(nrow(ell), 37L)
  expect_equal(ell[1, ], ell[37, ], tolerance = 1e-6)
})

## ---- as_xy helper ----

test_that("as_xy handles various inputs", {
  m <- cbind(1:3, 4:6)
  expect_equal(nrow(tissot:::as_xy(m)), 3L)
  expect_equal(nrow(tissot:::as_xy(data.frame(a = 1, b = 2))), 1L)
  expect_equal(nrow(tissot:::as_xy(c(1, 2))), 1L)
  expect_equal(nrow(tissot:::as_xy(list(x = 1:2, y = 3:4))), 2L)
  expect_equal(nrow(tissot:::as_xy(list(lon = 147, lat = -42))), 1L)
})

## ---- resolve_gpar ----

test_that("resolve_gpar works for TRUE, FALSE, and list", {
  defs <- list(col = "red", lwd = 2)
  expect_equal(tissot:::resolve_gpar(TRUE, defs), defs)
  expect_null(tissot:::resolve_gpar(FALSE, defs))
  merged <- tissot:::resolve_gpar(list(col = "blue"), defs)
  expect_equal(merged$col, "blue")
  expect_equal(merged$lwd, 2)
})

## ---- multiple projection types ----

test_that("multiple projection families work", {
  projs <- c("+proj=robin", "+proj=moll", "+proj=merc",
             "+proj=laea +lat_0=-90", "+proj=stere +lat_0=90",
             "+proj=utm +zone=55 +south")
  for (p in projs) {
    r <- tissot(cbind(147, -42), p)
    expect_true(all(is.finite(r$scale.area)),
                info = paste("Failed for", p))
  }
})
