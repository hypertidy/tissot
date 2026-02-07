test_that("conformal projection has near-zero angular deformation", {
  ## Stereographic is conformal
  r <- tissot(c(0, 45), "+proj=stere +lat_0=90")
  expect_lt(r$angle_deformation, 0.01)
  ## scale.a and scale.b should be equal for conformal
  expect_equal(r$scale.a, r$scale.b, tolerance = 1e-4)
})

test_that("equal-area projection has unit areal scale", {
  ## Lambert Azimuthal Equal Area
  r <- tissot(c(0, 45), "+proj=laea +lat_0=45")
  expect_equal(r$scale.area, 1.0, tolerance = 1e-4)
})

test_that("Plate Carree at equator has unit scale factors", {
  r <- tissot(c(0, 0), "+proj=eqc")
  expect_equal(r$scale.h, 1.0, tolerance = 0.01)
  expect_equal(r$scale.k, 1.0, tolerance = 0.01)
})

test_that("vectorized call matches single-point calls", {
  pts <- cbind(c(0, 30, -60), c(0, 45, -30))
  proj <- "+proj=robin"

  r_vec <- tissot(pts, proj)

  for (i in seq_len(nrow(pts))) {
    r_single <- tissot(pts[i, ], proj)
    expect_equal(r_vec$scale.a[i], r_single$scale.a, tolerance = 1e-10)
    expect_equal(r_vec$scale.b[i], r_single$scale.b, tolerance = 1e-10)
    expect_equal(r_vec$scale.area[i], r_single$scale.area, tolerance = 1e-10)
    expect_equal(r_vec$angle_deformation[i], r_single$angle_deformation,
                 tolerance = 1e-10)
  }
})

test_that("matrix input works the same as data.frame input", {
  m <- cbind(c(10, 20, 30), c(-10, 0, 10))
  proj <- "+proj=moll"

  r_mat <- tissot(m, proj)
  r_df <- tissot(as.data.frame(m), proj)

  expect_equal(r_mat$scale.a, r_df$scale.a)
  expect_equal(r_mat$scale.area, r_df$scale.area)
  expect_equal(r_mat$convergence, r_df$convergence)
})

test_that("tissot returns a tissot_tbl with correct attributes", {
  r <- tissot(c(0, 45), "+proj=robin")
  expect_s3_class(r, "tissot_tbl")
  expect_s3_class(r, "tbl_df")
  expect_equal(attr(r, "source"), "EPSG:4326")
  expect_equal(attr(r, "target"), "+proj=robin")
})

test_that("tissot_tbl has all expected columns", {
  r <- tissot(c(0, 45), "+proj=robin")
  expected <- c("x", "y", "dx_dlam", "dy_dlam", "dx_dphi", "dy_dphi",
                "scale.h", "scale.k", "scale.omega", "scale.a", "scale.b",
                "scale.area", "angle_deformation", "convergence")
  expect_true(all(expected %in% names(r)))
})

test_that("indicatrix accepts a tissot_tbl", {
  r <- tissot(cbind(c(0, 30), c(0, 45)), "+proj=robin")
  ii <- indicatrix(r)
  expect_s3_class(ii, "indicatrix_list")
  expect_length(ii, 2L)
  expect_s3_class(ii[[1]], "indicatrix")
  expect_equal(attr(ii, "target"), "+proj=robin")
})

test_that("indicatrix accepts raw coords with target", {
  ii <- indicatrix(c(0, 45), "+proj=stere +lat_0=90")
  expect_s3_class(ii, "indicatrix_list")
  expect_length(ii, 1L)
})

test_that("indicatrix errors when target is missing", {
  expect_error(indicatrix(c(0, 45)), "target is required")
})

test_that("ti_ellipse returns a closed polygon", {
  ii <- indicatrix(c(0, 45), "+proj=robin")
  ell <- ti_ellipse(ii[[1]], n = 36)
  ## First row should equal last row (closed polygon)
  expect_equal(ell[1L, ], ell[nrow(ell), ], tolerance = 1e-10)
  expect_equal(nrow(ell), 37L)
})

test_that("indicatrix_list subsetting preserves class", {
  r <- tissot(cbind(c(0, 30, 60), c(0, 30, 60)), "+proj=robin")
  ii <- indicatrix(r)
  sub <- ii[1:2]
  expect_s3_class(sub, "indicatrix_list")
  expect_length(sub, 2L)
  expect_equal(attr(sub, "target"), attr(ii, "target"))
})

test_that("multiple conformal projections all show near-zero deformation", {
  projs <- c(
    "+proj=stere +lat_0=90",
    "+proj=merc",
    "+proj=tmerc +lon_0=147"
  )
  for (p in projs) {
    r <- tissot(c(0, 30), p)
    expect_lt(r$angle_deformation, 0.05,
              label = sprintf("deformation for %s", p))
  }
})

test_that("multiple equal-area projections all have unit area", {
  projs <- c(
    "+proj=laea +lat_0=0",
    "+proj=moll",
    "+proj=sinu"
  )
  for (p in projs) {
    r <- tissot(c(0, 30), p)
    expect_equal(r$scale.area, 1.0, tolerance = 5e-3,
                 label = sprintf("area for %s", p))
  }
})

test_that("polar points don't produce NaN", {
  ## South pole in LAEA centered there
  r <- tissot(c(0, -89.99), "+proj=laea +lat_0=-90")
  expect_false(any(is.nan(r$scale.a)))
  expect_false(any(is.nan(r$scale.area)))
})

test_that("scale.omega equals angle_deformation", {
  r <- tissot(cbind(seq(-60, 60, by = 30), 0), "+proj=robin")
  expect_equal(r$scale.omega, r$angle_deformation)
})

test_that("print and summary methods work", {
  r <- tissot(cbind(c(0, 30), c(0, 45)), "+proj=robin")
  expect_output(print(r), "Tissot indicatrix: 2 points")
  expect_output(summary(r), "Target CRS:.*robin")

  ii <- indicatrix(r)
  expect_output(print(ii), "Indicatrix list: 2 ellipses")
})

test_that("custom source CRS works", {
  r <- tissot(c(0, 45), "+proj=robin", source = "EPSG:4267")
  expect_equal(nrow(r), 1L)
  expect_true(is.finite(r$scale.area))
  expect_equal(attr(r, "source"), "EPSG:4267")
})

test_that("as_xy handles various inputs", {
  m <- cbind(1:3, 4:6)
  expect_equal(nrow(tissot:::as_xy(m)), 3L)
  expect_equal(nrow(tissot:::as_xy(data.frame(a = 1, b = 2))), 1L)
  expect_equal(nrow(tissot:::as_xy(c(1, 2))), 1L)
  expect_equal(nrow(tissot:::as_xy(list(x = 1:2, y = 3:4))), 2L)
  expect_equal(nrow(tissot:::as_xy(list(lon = 147, lat = -42))), 1L)
})

test_that("resolve_gpar works for TRUE, FALSE, and list", {
  defs <- list(col = "red", lwd = 2)
  expect_equal(tissot:::resolve_gpar(TRUE, defs), defs)
  expect_null(tissot:::resolve_gpar(FALSE, defs))
  merged <- tissot:::resolve_gpar(list(col = "blue"), defs)
  expect_equal(merged$col, "blue")
  expect_equal(merged$lwd, 2)
})

test_that("tissot_map returns the world data invisibly", {
  skip_if_not(interactive(), "tissot_map needs a graphics device")
  r <- tissot(c(0, 0), "+proj=robin")
  ii <- indicatrix(r)
  plot(ii, add = FALSE)
  w <- tissot_map()
  expect_true(is.matrix(w))
  expect_equal(ncol(w), 2L)
})
