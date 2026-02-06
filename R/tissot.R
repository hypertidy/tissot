## Attempt to formalize Bill Huber's Tissot indicatrix calculations
## https://gis.stackexchange.com/a/5075/482

#' Compute Tissot indicatrix properties
#'
#' Compute the Tissot indicatrix at given longitude/latitude locations for a
#' map projection. Returns scale factors, angular deformation, convergence,
#' and related distortion properties.
#'
#' The Jacobian of the projection is computed via finite differences, projecting
#' all base and offset points in a single batched call to
#' [gdalraster::transform_xy()]. All subsequent calculations (SVD, distortion
#' metrics) are fully vectorized.
#'
#' @param x longitude values (or a two-column matrix/data.frame of lon, lat)
#' @param y latitude values (ignored if x is a matrix)
#' @param proj.in input CRS (default "EPSG:4326")
#' @param proj.out target projection CRS
#' @param A semi-major axis of the ellipsoid (default WGS84)
#' @param f.inv inverse flattening (default WGS84)
#' @param dx finite difference step in degrees (default 1e-4)
#' @return A tibble with columns: x (lon), y (lat), dx_dlam, dy_dlam,
#'   dx_dphi, dy_dphi, scale.h, scale.k, scale.omega, scale.a,
#'   scale.b, scale.area, angle_deformation, convergence
#' @export
#' @examples
#' tissot(0, 45, proj.out = "+proj=robin")
#' tissot(seq(-180, 180, by = 30), 0, proj.out = "+proj=robin")
tissot <- function(x, y = NULL,
                   proj.in = "EPSG:4326",
                   proj.out,
                   A = 6378137,
                   f.inv = 298.257223563,
                   dx = 1e-4) {

  ## Handle matrix/data.frame input
  if (is.matrix(x) || is.data.frame(x)) {
    y <- x[, 2L]
    x <- x[, 1L]
  }
  stopifnot(length(x) == length(y))

  n <- length(x)
  lon <- as.double(x)
  lat <- as.double(y)

  ## Ellipsoid derived quantities
  e.sq <- (2 - 1 / f.inv) / f.inv
  phi <- lat * pi / 180

  ## Radii of curvature (Snyder 4-18 to 4-21)
  sin_phi <- sin(phi)
  cos_phi <- cos(phi)
  sin2 <- sin_phi^2
  N_radius <- A / sqrt(1 - e.sq * sin2)
  R_meridian <- A * (1 - e.sq) / (1 - e.sq * sin2)^1.5

  ## Metric scale: degrees to metres on the ellipsoid
  h_phi <- R_meridian * (pi / 180)
  h_lam <- N_radius * cos_phi * (pi / 180)

  ## Build the 3N-point matrix: [base; lam+dx; phi+dx]
  pts <- matrix(NA_real_, nrow = 3L * n, ncol = 2L)
  pts[seq_len(n), ] <- cbind(lon, lat)
  pts[n + seq_len(n), ] <- cbind(lon + dx, lat)
  pts[2L * n + seq_len(n), ] <- cbind(lon, lat + dx)

  ## Single batched projection call
  proj <- gdalraster::transform_xy(pts, proj.in, proj.out)

  ## Extract projected coordinates
  idx_base <- seq_len(n)
  idx_lam <- n + seq_len(n)
  idx_phi <- 2L * n + seq_len(n)

  X0 <- proj[idx_base, 1L]
  Y0 <- proj[idx_base, 2L]
  X_lam <- proj[idx_lam, 1L]
  Y_lam <- proj[idx_lam, 2L]
  X_phi <- proj[idx_phi, 1L]
  Y_phi <- proj[idx_phi, 2L]

  ## Jacobian in projected coords per degree
  dX_dlam_deg <- (X_lam - X0) / dx
  dY_dlam_deg <- (Y_lam - Y0) / dx
  dX_dphi_deg <- (X_phi - X0) / dx
  dY_dphi_deg <- (Y_phi - Y0) / dx

  ## Jacobian normalised to dimensionless scale factors
  ## Dividing by h_lam/h_phi converts from projected-per-degree
  ## to projected-per-metre-on-the-ellipsoid
  dX_dlam <- dX_dlam_deg / h_lam
  dY_dlam <- dY_dlam_deg / h_lam
  dX_dphi <- dX_dphi_deg / h_phi
  dY_dphi <- dY_dphi_deg / h_phi

  ## Tissot parameters from the Jacobian (vectorized)
  ## Components of J^T J (2x2 symmetric matrix per point)
  a11 <- dX_dlam^2 + dY_dlam^2
  a12 <- dX_dlam * dX_dphi + dY_dlam * dY_dphi
  a22 <- dX_dphi^2 + dY_dphi^2

  ## Eigenvalues of J^T J via quadratic formula
  tr <- a11 + a22
  det <- a11 * a22 - a12^2
  disc <- sqrt(pmax(tr^2 - 4 * det, 0))

  lam1 <- (tr + disc) / 2
  lam2 <- (tr - disc) / 2

  ## Scale factors: singular values of J
  scale.a <- sqrt(pmax(lam1, 0))
  scale.b <- sqrt(pmax(lam2, 0))

  ## Meridional scale (h) and parallel scale (k) from Snyder
  scale.h <- sqrt(a22)
  scale.k <- sqrt(a11)

  ## Areal scale
  scale.area <- scale.a * scale.b

  ## Maximum angular deformation (omega): Snyder eq 4-1b
  sin_half_omega <- ifelse(scale.a + scale.b > 0,
                           (scale.a - scale.b) / (scale.a + scale.b),
                           0)
  angle_deformation <- 2 * asin(pmin(sin_half_omega, 1)) * 180 / pi

  ## Convergence: angle of the projected meridian from north
  convergence <- atan2(dX_dphi_deg, dY_dphi_deg) * 180 / pi

  tibble::tibble(
    x = lon,
    y = lat,
    dx_dlam = dX_dlam,
    dy_dlam = dY_dlam,
    dx_dphi = dX_dphi,
    dy_dphi = dY_dphi,
    scale.h = scale.h,
    scale.k = scale.k,
    scale.omega = angle_deformation,
    scale.a = scale.a,
    scale.b = scale.b,
    scale.area = scale.area,
    angle_deformation = angle_deformation,
    convergence = convergence
  )
}


#' Compute Tissot indicatrix ellipse geometry
#'
#' Generates indicatrix objects for plotting at each point.
#' Returns a list of indicatrix objects suitable for [plot.indicatrix()].
#'
#' @param x longitude values (or a two-column matrix)
#' @param y latitude values
#' @param proj.out target projection CRS
#' @param proj.in input CRS (default "EPSG:4326")
#' @param ... passed to [tissot()]
#' @return A list of indicatrix objects, each containing ellipse geometry
#'   and distortion properties
#' @export
indicatrix <- function(x, y = NULL, proj.out, proj.in = "EPSG:4326", ...) {
  if (is.matrix(x) || is.data.frame(x)) {
    y <- x[, 2L]
    x <- x[, 1L]
  }

  tis <- tissot(x, y, proj.in = proj.in, proj.out = proj.out, ...)

  ## Project the center points
  centers <- gdalraster::transform_xy(
    cbind(tis$x, tis$y), proj.in, proj.out
  )

  lapply(seq_len(nrow(tis)), function(i) {
    structure(list(
      center = centers[i, ],
      A = matrix(c(tis$dx_dlam[i], tis$dy_dlam[i],
                    tis$dx_dphi[i], tis$dy_dphi[i]), 2L, 2L),
      scale.h = tis$scale.h[i],
      scale.k = tis$scale.k[i],
      scale.a = tis$scale.a[i],
      scale.b = tis$scale.b[i],
      scale.area = tis$scale.area[i],
      angle_deformation = tis$angle_deformation[i],
      convergence = tis$convergence[i],
      lon = tis$x[i],
      lat = tis$y[i]
    ), class = "indicatrix")
  })
}


#' Generate ellipse coordinates for an indicatrix
#'
#' @param x an indicatrix object
#' @param scale scaling factor for the ellipse size in projected units
#' @param n number of points on the ellipse
#' @param ... ignored
#' @return A two-column matrix of projected coordinates tracing the ellipse
#' @export
ti_ellipse <- function(x, scale = 1e5, n = 72, ...) {
  stopifnot(inherits(x, "indicatrix"))
  theta <- seq(0, 2 * pi, length.out = n + 1L)
  circle <- cbind(cos(theta), sin(theta))
  ellipse <- circle %*% t(x$A) * scale
  ellipse[, 1L] <- ellipse[, 1L] + x$center[1L]
  ellipse[, 2L] <- ellipse[, 2L] + x$center[2L]
  ellipse
}


#' Plot an indicatrix
#'
#' @param x an indicatrix object
#' @param scale scaling factor for the ellipse size
#' @param n number of points on the ellipse
#' @param col fill colour for the ellipse
#' @param border border colour
#' @param add logical; add to existing plot?
#' @param ... passed to [polygon()]
#' @export
plot.indicatrix <- function(x, scale = 1e5, n = 72,
                            col = "#FF990055", border = "black",
                            add = TRUE, ...) {
  ell <- ti_ellipse(x, scale = scale, n = n)
  if (!add) {
    plot(ell, type = "n", asp = 1, xlab = "", ylab = "")
  }
  polygon(ell[, 1L], ell[, 2L], col = col, border = border, ...)
  invisible(x)
}
