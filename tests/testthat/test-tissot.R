test_that("conformal projection has near-zero angular deformation", {
  ## Stereographic is conformal
  r <- tissot(0, 45, proj.out = "+proj=stere +lat_0=90")
  expect_lt(r$angle_deformation, 0.01)
  ## scale.a and scale.b should be equal for conformal

  expect_equal(r$scale.a, r$scale.b, tolerance = 1e-4)
})

test_that("equal-area projection has unit areal scale", {
  ## Lambert Azimuthal Equal Area
  r <- tissot(0, 45, proj.out = "+proj=laea +lat_0=45")
  expect_equal(r$scale.area, 1.0, tolerance = 1e-4)
})

test_that("Plate Carree at equator has unit scale factors",
          {
            r <- tissot(0, 0, proj.out = "+proj=eqc")
            expect_equal(r$scale.h, 1.0, tolerance = 0.01)
            expect_equal(r$scale.k, 1.0, tolerance = 0.01)
          })

test_that("vectorized call matches single-point calls", {
  lons <- c(0, 30, -60)
  lats <- c(0, 45, -30)
  proj <- "+proj=robin"

  r_vec <- tissot(lons, lats, proj.out = proj)

  for (i in seq_along(lons)) {
    r_single <- tissot(lons[i], lats[i], proj.out = proj)
    expect_equal(r_vec$scale.a[i], r_single$scale.a, tolerance = 1e-10)
    expect_equal(r_vec$scale.b[i], r_single$scale.b, tolerance = 1e-10)
    expect_equal(r_vec$scale.area[i], r_single$scale.area, tolerance = 1e-10)
    expect_equal(r_vec$angle_deformation[i], r_single$angle_deformation,
                 tolerance = 1e-10)
  }
})

test_that("matrix input works the same as separate vectors", {
  lons <- c(10, 20, 30)
  lats <- c(-10, 0, 10)
  proj <- "+proj=moll"

  r_xy <- tissot(lons, lats, proj.out = proj)
  r_mat <- tissot(cbind(lons, lats), proj.out = proj)

  expect_equal(r_xy$scale.a, r_mat$scale.a)
  expect_equal(r_xy$scale.area, r_mat$scale.area)
  expect_equal(r_xy$convergence, r_mat$convergence)
})

test_that("tissot returns a tissot_tbl with correct attributes", {
  r <- tissot(0, 45, proj.out = "+proj=robin")
  expect_s3_class(r, "tissot_tbl")
  expect_s3_class(r, "tbl_df")
  expect_equal(attr(r, "proj.in"), "EPSG:4326")
  expect_equal(attr(r, "proj.out"), "+proj=robin")
})

test_that("tissot_tbl has all expected columns", {
  r <- tissot(0, 45, proj.out = "+proj=robin")
  expected <- c("x", "y", "dx_dlam", "dy_dlam", "dx_dphi", "dy_dphi",
                "scale.h", "scale.k", "scale.omega", "scale.a", "scale.b",
                "scale.area", "angle_deformation", "convergence")
  expect_true(all(expected %in% names(r)))
})

test_that("indicatrix accepts a tissot_tbl (item 3)", {
  r <- tissot(c(0, 30), c(0, 45), proj.out = "+proj=robin")
  ii <- indicatrix(r)
  expect_s3_class(ii, "indicatrix_list")
  expect_length(ii, 2L)
  expect_s3_class(ii[[1]], "indicatrix")
  expect_equal(attr(ii, "proj.out"), "+proj=robin")
})

test_that("indicatrix accepts raw lon/lat with proj.out", {
  ii <- indicatrix(0, 45, proj.out = "+proj=stere +lat_0=90")
  expect_s3_class(ii, "indicatrix_list")
  expect_length(ii, 1L)
})

test_that("indicatrix errors when proj.out is missing", {
  expect_error(indicatrix(0, 45), "proj.out is required")
})

test_that("ti_ellipse returns a closed polygon", {
  ii <- indicatrix(0, 45, proj.out = "+proj=robin")
  ell <- ti_ellipse(ii[[1]], n = 36)
  ## First row should equal last row (closed polygon)
  expect_equal(ell[1L, ], ell[nrow(ell), ], tolerance = 1e-10)
  expect_equal(nrow(ell), 37L)
})

test_that("indicatrix_list subsetting preserves class", {
  r <- tissot(c(0, 30, 60), c(0, 30, 60), proj.out = "+proj=robin")
  ii <- indicatrix(r)
  sub <- ii[1:2]
  expect_s3_class(sub, "indicatrix_list")
  expect_length(sub, 2L)
  expect_equal(attr(sub, "proj.out"), attr(ii, "proj.out"))
})

test_that("multiple conformal projections all show near-zero deformation", {
  projs <- c(
    "+proj=stere +lat_0=90",
    "+proj=merc",
    "+proj=tmerc +lon_0=147"
  )
  for (p in projs) {
    r <- tissot(0, 30, proj.out = p)
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
    r <- tissot(0, 30, proj.out = p)
    expect_equal(r$scale.area, 1.0, tolerance = 5e-3,
                 label = sprintf("area for %s", p))
  }
})

test_that("polar points don't produce NaN", {
  ## South pole in LAEA centered there
  r <- tissot(0, -89.99, proj.out = "+proj=laea +lat_0=-90")
  expect_false(any(is.nan(r$scale.a)))
  expect_false(any(is.nan(r$scale.area)))
})

test_that("scale.omega equals angle_deformation", {
  r <- tissot(cbind(seq(-60, 60, by = 30), 0), proj.out = "+proj=robin")
  expect_equal(r$scale.omega, r$angle_deformation)
})

test_that("print and summary methods work", {
  r <- tissot(c(0, 30), c(0, 45), proj.out = "+proj=robin")
  expect_output(print(r), "Tissot indicatrix: 2 points")
  expect_output(summary(r), "Target CRS:.*robin")

  ii <- indicatrix(r)
  expect_output(print(ii), "Indicatrix list: 2 ellipses")
})

test_that("tissot_map returns the world data invisibly", {
  skip_if_not(interactive(), "tissot_map needs a graphics device")
  r <- tissot(0, 0, proj.out = "+proj=robin")
  ii <- indicatrix(r)
  plot(ii, add = FALSE)
  w <- tissot_map()
  expect_true(is.matrix(w))
  expect_equal(ncol(w), 2L)
})
