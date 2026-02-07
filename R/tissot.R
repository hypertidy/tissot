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
#' @param x input coordinates — any xy-ish object: a two-column matrix,
#'   data.frame, tibble, list with `x`/`y` or `lon`/`lat` components,
#'   or a length-2 numeric vector for a single point
#' @param target target projection CRS string (required)
#' @param ... ignored
#' @param source source CRS (default `"EPSG:4326"`)
#' @param A semi-major axis of the ellipsoid (default WGS84)
#' @param f.inv inverse flattening (default WGS84)
#' @param dx finite difference step in degrees (default 1e-4)
#' @return A `tissot_tbl` tibble with columns: x (lon), y (lat), dx_dlam,
#'   dy_dlam, dx_dphi, dy_dphi, scale.h, scale.k, scale.omega, scale.a,
#'   scale.b, scale.area, angle_deformation, convergence. The `source` and
#'   `target` CRS strings are stored as attributes.
#' @seealso [indicatrix()], [tissot_map()]
#' @export
#' @examples
#' tissot(c(0, 45), "+proj=robin")
#' tissot(cbind(seq(-180, 180, by = 30), 0), "+proj=robin")
tissot <- function(x, target, ...,
                   source = "EPSG:4326",
                   A = 6378137,
                   f.inv = 298.257223563,
                   dx = 1e-4) {

  xy <- as_xy(x)
  n <- nrow(xy)
  lon <- xy[, 1L]
  lat <- xy[, 2L]

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
  proj <- gdalraster::transform_xy(pts, source, target)

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

  out <- tibble::tibble(
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

  ## Attach projection metadata and subclass
  attr(out, "source") <- source
  attr(out, "target") <- target
  class(out) <- c("tissot_tbl", class(out))
  out
}


#' @export
print.tissot_tbl <- function(x, ...) {
  target <- attr(x, "target") %||% "unknown"
  cat(sprintf("Tissot indicatrix: %d point%s, %s\n",
              nrow(x), if (nrow(x) != 1L) "s" else "", target))
  NextMethod()
}

#' @export
summary.tissot_tbl <- function(object, ...) {
  src <- attr(object, "source") %||% "unknown"
  tgt <- attr(object, "target") %||% "unknown"
  cat(sprintf("Tissot indicatrix: %d point%s\n", nrow(object),
              if (nrow(object) != 1L) "s" else ""))
  cat(sprintf("  Source CRS: %s\n", src))
  cat(sprintf("  Target CRS: %s\n", tgt))
  cat(sprintf("  Areal scale:  min=%.4f  max=%.4f  mean=%.4f\n",
              min(object$scale.area, na.rm = TRUE),
              max(object$scale.area, na.rm = TRUE),
              mean(object$scale.area, na.rm = TRUE)))
  cat(sprintf("  Angular def:  min=%.4f  max=%.4f  mean=%.4f deg\n",
              min(object$angle_deformation, na.rm = TRUE),
              max(object$angle_deformation, na.rm = TRUE),
              mean(object$angle_deformation, na.rm = TRUE)))
  cat(sprintf("  Scale h:      min=%.4f  max=%.4f  (meridional)\n",
              min(object$scale.h, na.rm = TRUE),
              max(object$scale.h, na.rm = TRUE)))
  cat(sprintf("  Scale k:      min=%.4f  max=%.4f  (parallel)\n",
              min(object$scale.k, na.rm = TRUE),
              max(object$scale.k, na.rm = TRUE)))
  invisible(object)
}


#' Compute Tissot indicatrix ellipse geometry
#'
#' Generates indicatrix objects for plotting at each point.
#' Returns an `indicatrix_list` containing individual indicatrix objects.
#'
#' `indicatrix()` accepts either:
#' - A `tissot_tbl` object (from [tissot()]) — projection is extracted
#'   from attributes, `target` is optional
#' - Any xy-ish input with an explicit `target`
#'
#' @param x a `tissot_tbl`, or any xy-ish input (see [tissot()])
#' @param target target projection CRS (extracted from `tissot_tbl`
#'   attributes if `x` is one; required otherwise)
#' @param ... passed to [tissot()]
#' @param source source CRS (default `"EPSG:4326"`)
#' @return An `indicatrix_list` object (a list of `indicatrix` objects with
#'   `source` and `target` stored as attributes)
#' @seealso [tissot()], [plot.indicatrix()], [plot.indicatrix_list()],
#'   [ti_ellipse()]
#' @export
#' @examples
#' ## From a tissot_tbl
#' r <- tissot(cbind(seq(-150, 150, by = 30), 0), "+proj=robin")
#' ii <- indicatrix(r)
#' plot(ii, scale = 5e5)
#'
#' ## From raw coordinates
#' ii2 <- indicatrix(c(0, 45), "+proj=stere +lat_0=90")
#' plot(ii2)
indicatrix <- function(x, target = NULL, ..., source = "EPSG:4326") {

  if (inherits(x, "tissot_tbl")) {
    ## Extract projection from the tissot_tbl attributes
    tis <- x
    source <- attr(tis, "source") %||% source
    target <- target %||% attr(tis, "target")
    if (is.null(target)) {
      stop("target not found in tissot_tbl attributes and not supplied")
    }
  } else {
    ## Raw coordinate path
    if (is.null(target)) {
      stop("target is required when x is not a tissot_tbl")
    }
    tis <- tissot(x, target, source = source, ...)
  }

  ## Project the center points
  centers <- gdalraster::transform_xy(
    cbind(tis$x, tis$y), source, target
  )

  out <- lapply(seq_len(nrow(tis)), function(i) {
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

  structure(out,
            class = "indicatrix_list",
            source = source,
            target = target)
}


#' Generate ellipse coordinates for an indicatrix
#'
#' @param x an indicatrix object
#' @param scale scaling factor for the ellipse size in projected units
#' @param n number of points on the ellipse
#' @param ... ignored
#' @return A two-column matrix of projected coordinates tracing the ellipse
#' @seealso [indicatrix()], [plot.indicatrix()]
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


#' Generate reference unit circle coordinates
#'
#' The unit circle shows what a circle of unit scale factor looks like —
#' comparing the indicatrix ellipse to this circle shows the distortion.
#'
#' @param x an indicatrix object
#' @param scale scaling factor (should match the indicatrix plot scale)
#' @param n number of points on the circle
#' @return A two-column matrix of projected coordinates
#' @keywords internal
#' @noRd
ti_circle <- function(x, scale = 1e5, n = 72) {
  theta <- seq(0, 2 * pi, length.out = n + 1L)
  circle <- cbind(cos(theta), sin(theta)) * scale
  circle[, 1L] <- circle[, 1L] + x$center[1L]
  circle[, 2L] <- circle[, 2L] + x$center[2L]
  circle
}


#' Generate axis line coordinates for an indicatrix
#'
#' Computes the endpoints for the lambda (parallel) and phi (meridian)
#' direction lines through the indicatrix center.
#'
#' @param x an indicatrix object
#' @param scale scaling factor
#' @return A list with `lambda` and `phi` components, each a 2x2 matrix
#'   of endpoint coordinates
#' @keywords internal
#' @noRd
ti_axes <- function(x, scale = 1e5) {
  ## Lambda direction (parallel scale, columns 1 of A)
  lam_dir <- x$A[, 1L]
  lam_dir <- lam_dir / sqrt(sum(lam_dir^2))
  lam_len <- x$scale.k * scale

  lam <- rbind(
    x$center - lam_dir * lam_len,
    x$center + lam_dir * lam_len
  )

  ## Phi direction (meridional scale, columns 2 of A)
  phi_dir <- x$A[, 2L]
  phi_dir <- phi_dir / sqrt(sum(phi_dir^2))
  phi_len <- x$scale.h * scale

  phi <- rbind(
    x$center - phi_dir * phi_len,
    x$center + phi_dir * phi_len
  )

  list(lambda = lam, phi = phi)
}


#' Plot an indicatrix
#'
#' Draws a single Tissot indicatrix ellipse on the current plot. The ellipse
#' shows the distortion of a unit circle under the map projection. Optional
#' overlays include a reference unit circle, and lambda/phi direction axes.
#'
#' `show.circle` and `show.axes` accept `TRUE` (use defaults), `FALSE`
#' (hide), or a named list of graphical parameters to override defaults.
#' For example, `show.circle = list(border = "blue", lty = 3)`.
#'
#' @param x an `indicatrix` object (from [indicatrix()])
#' @param scale scaling factor for the ellipse size in projected units
#' @param n number of points on the ellipse
#' @param col fill colour for the ellipse
#' @param border border colour
#' @param add logical; add to existing plot?
#' @param show.axes `TRUE`, `FALSE`, or a named list of graphical parameters
#'   for the direction lines. Defaults: `list(col.lambda = "red",
#'   col.phi = "blue", lwd = 1.5)`.
#' @param show.circle `TRUE`, `FALSE`, or a named list of graphical parameters
#'   for the reference circle. Defaults: `list(col = adjustcolor("white",
#'   alpha.f = 0.6), border = "grey70", lwd = 2.5, lty = 2)`.
#' @param ... passed to [graphics::polygon()]
#' @seealso [indicatrix()], [plot.indicatrix_list()], [ti_ellipse()]
#' @export
plot.indicatrix <- function(x, scale = 1e5, n = 72,
                            col = "#FF990055", border = "black",
                            add = TRUE,
                            show.axes = TRUE,
                            show.circle = TRUE, ...) {
  ell <- ti_ellipse(x, scale = scale, n = n)

  if (!add) {
    plot(ell, type = "n", asp = 1, xlab = "", ylab = "")
  }

  ## Reference unit circle (behind the ellipse)
  circle_gp <- resolve_gpar(show.circle, list(
    col = grDevices::adjustcolor("white", alpha.f = 1),
    border = "grey70", lwd = 4.5,
    lty = 2
  ))
  if (!is.null(circle_gp)) {
    circ <- ti_circle(x, scale = scale, n = n)
    graphics::polygon(circ[, 1L], circ[, 2L],
                      col = circle_gp$col,
                      border = circle_gp$border,
                      lty = circle_gp$lty, lwd = circle_gp$lwd)
  }

  ## Filled ellipse
  graphics::polygon(ell[, 1L], ell[, 2L], col = col, border = border, ...)

  ## Direction axes
  axes_gp <- resolve_gpar(show.axes, list(
    col.lambda = "red",
    col.phi = "blue",
    lwd = 1.5
  ))
  if (!is.null(axes_gp)) {
    ax <- ti_axes(x, scale = scale)
    graphics::lines(ax$lambda[, 1L], ax$lambda[, 2L],
                    col = axes_gp$col.lambda, lwd = axes_gp$lwd)
    graphics::lines(ax$phi[, 1L], ax$phi[, 2L],
                    col = axes_gp$col.phi, lwd = axes_gp$lwd)
  }

  invisible(x)
}


#' Plot a list of indicatrixes
#'
#' Draws all indicatrixes in an `indicatrix_list`, optionally creating a
#' new plot or adding to an existing one. Can colour-code the fill by a
#' distortion metric.
#'
#' @param x an `indicatrix_list` (from [indicatrix()])
#' @param scale scaling factor for ellipse size in projected units
#' @param n number of points per ellipse
#' @param col fill colour. If a single colour, used for all ellipses. If
#'   `NULL` and `fill.by` is set, colours are generated automatically.
#' @param border border colour
#' @param add logical; add to existing plot? If `FALSE`, creates a new
#'   plot sized to contain all ellipses.
#' @param show.axes `TRUE`, `FALSE`, or a named list of graphical parameters.
#'   See [plot.indicatrix()] for defaults.
#' @param show.circle `TRUE`, `FALSE`, or a named list of graphical parameters.
#'   Default `TRUE` for the list method (the circle-vs-ellipse comparison makes
#'   distortion visible at map scale). See [plot.indicatrix()] for defaults.
#' @param fill.by character; name of a distortion metric to colour-code
#'   the fill. One of `"scale.area"`, `"angle_deformation"`, `"scale.h"`,
#'   `"scale.k"`, `"scale.a"`, `"scale.b"`. Default `NULL` (uniform fill).
#' @param palette colour palette function (default [grDevices::hcl.colors()])
#' @param ncolors number of colours in the palette (default 64)
#' @param ... passed to [plot.indicatrix()]
#' @seealso [indicatrix()], [plot.indicatrix()], [tissot_map()]
#' @export
#' @examples
#' xy <- expand.grid(seq(-150, 150, by = 30), seq(-60, 60, by = 30))
#' r <- tissot(xy, "+proj=robin")
#' ii <- indicatrix(r)
#'
#' ## Uniform fill
#' plot(ii, scale = 6e5, add = FALSE)
#' tissot_map()
#'
#' ## Colour by areal distortion
#' plot(ii, scale = 6e5, add = FALSE, fill.by = "scale.area")
#' tissot_map()
plot.indicatrix_list <- function(x, scale = 1e5, n = 72,
                                 col = "#FF990055",
                                 border = "black",
                                 add = FALSE,
                                 show.axes = TRUE,
                                 show.circle = TRUE,
                                 fill.by = NULL,
                                 palette = NULL,
                                 ncolors = 64L, ...) {

  if (!add) {
    ## Compute bounding box from all ellipses
    all_pts <- do.call(rbind, lapply(x, ti_ellipse, scale = scale, n = n))
    xr <- range(all_pts[, 1L], na.rm = TRUE)
    yr <- range(all_pts[, 2L], na.rm = TRUE)
    plot(xr, yr, type = "n", asp = 1, xlab = "", ylab = "")

    ## Register the projection for tissot_map()
    target <- attr(x, "target")
    if (!is.null(target)) {
      options(tissot.last.plot.proj = target)
    }
  }

  ## Colour mapping
  cols <- rep_len(col, length(x))
  if (!is.null(fill.by)) {
    vals <- vapply(x, function(xi) xi[[fill.by]], numeric(1L))
    if (is.null(palette)) {
      palette <- function(n) grDevices::hcl.colors(n, palette = "YlOrRd",
                                                   rev = TRUE)
    }
    pal <- palette(ncolors)
    rng <- range(vals, na.rm = TRUE)
    if (diff(rng) > 0) {
      idx <- findInterval(vals, seq(rng[1L], rng[2L], length.out = ncolors + 1L),
                          all.inside = TRUE)
      cols <- grDevices::adjustcolor(pal[idx], alpha.f = 0.7)
    }
  }

  for (i in seq_along(x)) {
    plot.indicatrix(x[[i]], scale = scale, n = n,
                    col = cols[i], border = border, add = TRUE,
                    show.axes = show.axes, show.circle = show.circle, ...)
  }

  invisible(x)
}


#' Subset an indicatrix_list
#'
#' @param x an `indicatrix_list`
#' @param i indices
#' @return An `indicatrix_list` with the selected elements
#' @export
`[.indicatrix_list` <- function(x, i) {
  out <- unclass(x)[i]
  structure(out,
            class = "indicatrix_list",
            source = attr(x, "source"),
            target = attr(x, "target"))
}

#' @export
print.indicatrix_list <- function(x, ...) {
  target <- attr(x, "target") %||% "unknown"
  cat(sprintf("Indicatrix list: %d ellipse%s, %s\n",
              length(x), if (length(x) != 1L) "s" else "", target))
  invisible(x)
}

#' @export
length.indicatrix_list <- function(x) {
  length(unclass(x))
}
