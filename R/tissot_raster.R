#' Compute distortion surfaces on a regular grid
#'
#' Create a raster-like grid of Tissot distortion metrics in the target
#' projection. The grid is laid out in projected coordinates, back-projected
#' to longitude/latitude, and passed through [tissot()]. Returns a list of
#' matrices compatible with [image()] and easily converted to raster objects
#' (e.g. via [terra::rast()]).
#'
#' If `extent` is `NULL`, the function projects the four corners and midpoints
#' of the \eqn{[-180, 180] \times [-85, 85]}{[-180,180] x [-85,85]} lonlat
#' bounding box into the target CRS, and uses the resulting range (padded by
#' 5%) as the grid extent.
#'
#' @param target target projection CRS string
#' @param extent numeric length 4: `c(xmin, xmax, ymin, ymax)` in projected
#'   coordinates. If `NULL` (default), estimated from the global lonlat
#'   bounding box (clamped by `radius`).
#' @param radius maximum half-width of the grid extent in projected units
#'   (default `2e7`). The auto-detected extent is clamped so that no
#'   dimension exceeds `2 * radius`. Prevents blow-up for projections like
#'   stereographic or gnomonic that map the whole globe to huge coordinates.
#'   Ignored if `extent` is supplied.
#' @param nx integer; number of grid cells in x (default 100)
#' @param ny integer; number of grid cells in y. If `NULL`, set to maintain
#'   approximately square cells.
#' @param metrics character vector of column names from [tissot()] output to
#'   include as layers. Default includes areal scale, angular deformation,
#'   and meridional/parallel scale factors.
#' @param ... passed to [tissot()]
#' @param source source CRS (default `"EPSG:4326"`)
#' @return A `tissot_raster` list with components:
#'   \describe{
#'     \item{`x`}{numeric vector of x (easting) cell centers}
#'     \item{`y`}{numeric vector of y (northing) cell centers}
#'     \item{`z`}{named list of `ny x nx` matrices, one per metric}
#'     \item{`target`}{the target CRS string}
#'     \item{`extent`}{the `c(xmin, xmax, ymin, ymax)` used}
#'   }
#'   The object has an [image()] method for quick plotting.
#' @seealso [tissot()], [plot.tissot_raster()], [image.tissot_raster()]
#' @export
#' @examples
#' tr <- tissot_raster("+proj=robin", nx = 60)
#' image(tr)
#'
#' ## specific metric
#' image(tr, metric = "angle_deformation")
#'
#' ## with coastline overlay
#' image(tr)
#' tissot_map()
#'
#' ## polar stereo — use radius to cap extent
#' tp <- tissot_raster("+proj=stere +lat_0=-90", radius = 5e6, nx = 60)
#' image(tp)
#' tissot_map()
tissot_raster <- function(target, extent = NULL, radius = 2e7,
                          nx = 100L, ny = NULL,
                          metrics = c("scale.area", "angle_deformation",
                                      "scale.h", "scale.k"),
                          ..., source = "EPSG:4326") {

  ## Auto-detect extent from global bounding box, clamped by radius
  if (is.null(extent)) {
    extent <- .auto_extent(target, source, radius = radius)
  }
  stopifnot(length(extent) == 4L)
  xmin <- extent[1L]; xmax <- extent[2L]
  ymin <- extent[3L]; ymax <- extent[4L]

  ## Grid dimensions
  nx <- as.integer(nx)
  if (is.null(ny)) {
    ny <- as.integer(round(nx * (ymax - ymin) / (xmax - xmin)))
    ny <- max(ny, 1L)
  }
  ny <- as.integer(ny)

  ## Cell centers
  xc <- seq(xmin, xmax, length.out = nx)
  yc <- seq(ymin, ymax, length.out = ny)
  grid <- expand.grid(x = xc, y = yc)
  xy_proj <- as.matrix(grid)

  ## Back-project to lonlat
  xy_ll <- gdalraster::transform_xy(xy_proj, target, source)

  ## Flag valid back-projections (finite and within plausible lonlat range)
  ok <- is.finite(xy_ll[, 1L]) & is.finite(xy_ll[, 2L]) &
    abs(xy_ll[, 1L]) <= 360 & abs(xy_ll[, 2L]) <= 90

  ## Compute Tissot on valid cells
  z_list <- setNames(
    lapply(metrics, function(m) matrix(NA_real_, nrow = ny, ncol = nx)),
    metrics
  )

  if (any(ok)) {
    tis <- tissot(xy_ll[ok, , drop = FALSE], target, ..., source = source)

    ## Fill matrices — expand.grid fills x fastest, matching matrix[row, col]
    ## where row = y index, col = x index
    idx_ok <- which(ok)
    col_idx <- ((idx_ok - 1L) %% nx) + 1L   # x index
    row_idx <- ((idx_ok - 1L) %/% nx) + 1L   # y index

    for (m in metrics) {
      vals <- tis[[m]]
      if (!is.null(vals)) {
        z_list[[m]][cbind(row_idx, col_idx)] <- vals
      }
    }
  }

  ## Register projection for tissot_map()
  options(tissot.last.plot.proj = target)

  structure(list(
    x = xc,
    y = yc,
    z = z_list,
    target = target,
    extent = extent
  ), class = "tissot_raster")
}


#' @export
print.tissot_raster <- function(x, ...) {
  cat(sprintf("Tissot distortion raster: %d x %d, %s\n",
              length(x$x), length(x$y), x$target))
  cat(sprintf("  Extent: [%.0f, %.0f] x [%.0f, %.0f]\n",
              x$extent[1L], x$extent[2L], x$extent[3L], x$extent[4L]))
  cat(sprintf("  Metrics: %s\n", paste(names(x$z), collapse = ", ")))
  for (m in names(x$z)) {
    vals <- x$z[[m]]
    rng <- range(vals, na.rm = TRUE)
    if (all(is.finite(rng))) {
      cat(sprintf("    %s: [%.4f, %.4f]\n", m, rng[1L], rng[2L]))
    }
  }
  invisible(x)
}


#' Plot a Tissot distortion raster
#'
#' @param x a `tissot_raster` (from [tissot_raster()])
#' @param metric which metric to plot (default: first available)
#' @param col colour palette (default: `hcl.colors(64, "YlOrRd", rev = TRUE)`)
#' @param asp aspect ratio (default 1)
#' @param main plot title (default: metric name)
#' @param ... passed to [image()]
#' @rdname tissot_raster
#' @export
image.tissot_raster <- function(x, metric = NULL,
                                col = NULL,
                                asp = 1, main = NULL, ...) {
  metric <- metric %||% names(x$z)[1L]
  stopifnot(metric %in% names(x$z))
  if (is.null(col)) {
    col <- grDevices::hcl.colors(64L, palette = "YlOrRd", rev = TRUE)
  }
  if (is.null(main)) main <- metric

  ## Register projection for tissot_map()
  options(tissot.last.plot.proj = x$target)

  graphics::image(x$x, x$y, t(x$z[[metric]]),
                  col = col, asp = asp, main = main,
                  xlab = "", ylab = "", ...)
}


#' @rdname tissot_raster
#' @export
plot.tissot_raster <- function(x, ...) {
  image.tissot_raster(x, ...)
}


## ---- Internal: auto-detect projected extent ----

.auto_extent <- function(target, source = "EPSG:4326", radius = 2e7) {
  ## Sample points around the global bounding box
  lons <- c(-180, -90, 0, 90, 180, -180, -90, 0, 90, 180,
            -180, 0, 180, -180, 0, 180)
  lats <- c(-85, -85, -85, -85, -85, 0, 0, 0, 0, 0,
            85, 85, 85, -45, -45, 45)
  ## Add a denser ring to better capture curved projections
  ring_lon <- seq(-180, 180, by = 15)
  ring_lat_lo <- rep(-85, length(ring_lon))
  ring_lat_hi <- rep(85, length(ring_lon))
  ring_lat_mid <- rep(0, length(ring_lon))

  all_ll <- cbind(
    c(lons, ring_lon, ring_lon, ring_lon),
    c(lats, ring_lat_lo, ring_lat_hi, ring_lat_mid)
  )

  proj <- gdalraster::transform_xy(all_ll, source, target)
  ok <- is.finite(proj[, 1L]) & is.finite(proj[, 2L])

  if (sum(ok) < 2L) {
    stop("Cannot auto-detect extent for this projection. ",
         "Supply 'extent' manually.")
  }

  xr <- range(proj[ok, 1L])
  yr <- range(proj[ok, 2L])

  ## Pad 5%
  dx <- diff(xr) * 0.05
  dy <- diff(yr) * 0.05
  ext <- c(xr[1L] - dx, xr[2L] + dx, yr[1L] - dy, yr[2L] + dy)

  ## Clamp to radius around center
  cx <- (ext[1L] + ext[2L]) / 2
  cy <- (ext[3L] + ext[4L]) / 2
  ext[1L] <- max(ext[1L], cx - radius)
  ext[2L] <- min(ext[2L], cx + radius)
  ext[3L] <- max(ext[3L], cy - radius)
  ext[4L] <- min(ext[4L], cy + radius)

  ext
}
