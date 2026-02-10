#' Plot world coastline on a projected map
#'
#' `tissot_map()` draws the bundled [world] coastline, projected if a
#' projection is current. The projection is determined in this order:
#' 1. An explicit `target` argument
#' 2. The projection stored by the last [plot.indicatrix_list()] call
#'
#' `tissot_abline()` draws vertical and horizontal reference lines at a
#' given longitude/latitude in projected coordinates.
#'
#' @param target target CRS. If `NULL`, uses the last plot projection
#'   (from `plot.indicatrix_list()`) or draws in lon/lat.
#' @param add logical; add to existing plot (default `TRUE`) or create new
#' @param ... graphical parameters passed to [graphics::lines()] (if adding)
#'   or [graphics::plot()] (if creating new)
#' @return `tissot_map()` invisibly returns the (projected) world coastline matrix
#' @seealso [plot.indicatrix_list()], [tissot()]
#' @export
#' @examples
#' r <- tissot(cbind(seq(-150, 150, by = 30), 0), "+proj=robin")
#' ii <- indicatrix(r)
#' plot(ii, scale = 6e5, add = FALSE)
#' tissot_map()
tissot_map <- function(..., target = NULL, add = TRUE) {
  props <- list(...)

  if (is.null(target)) {
    target <- getOption("tissot.last.plot.proj")
  }

  w <- if (is.null(target)) {
    world
  } else {
    suppressWarnings(gdalraster::transform_xy(world, "OGC:CRS84", target))
  }

  if (is.null(props$col)) props$col <- grDevices::rgb(.7, .7, .7)
  if (add) {
    props$x <- w
    do.call(graphics::lines, props)
  } else {
    props$x <- w
    if (is.null(props$pch)) props$pch <- "."
    do.call(graphics::plot, props)
    if (!is.null(target)) {
      options(tissot.last.plot.proj = target)
    }
  }
  invisible(w)
}


#' @rdname tissot_map
#' @param x longitude values (or any xy-ish input; see [tissot()])
#' @param y latitude values (ignored if `x` is a matrix)
#' @param source source CRS for the coordinates (default `"EPSG:4326"`)
#' @return `tissot_abline()` is called for its side effect
#' @export
#' @importFrom graphics abline
tissot_abline <- function(x, y = NULL, ..., source = "EPSG:4326",
                          target = NULL) {
  xy <- as_xy(if (is.null(y)) x else cbind(x, y))

  if (is.null(target)) {
    target <- getOption("tissot.last.plot.proj")
  }
  if (!is.null(target)) {
    xy <- gdalraster::transform_xy(xy, source, target)
  }
  graphics::abline(v = xy[, 1L], h = xy[, 2L], ...)
}


#' Get or set the current plot projection
#'
#' @description
#'  **Deprecated.** Prefer passing `target` explicitly.
#'
#' These functions access a global option. Prefer passing `target`
#' explicitly to [tissot_map()] and [tissot_abline()], or let
#' [plot.indicatrix_list()] set the projection automatically.
#'
#' @return `tissot_get_proj()` returns the current projection string or `NULL`
#' @export
tissot_get_proj <- function() {
  getOption("tissot.last.plot.proj")
}

#' @rdname tissot_get_proj
#' @param target projection CRS string
#' @export
tissot_set_proj <- function(target) {
  options(tissot.last.plot.proj = target)
}
