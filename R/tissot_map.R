#' Plot world coastline on a projected map
#'
#' `tissot_map()` draws the bundled [world] coastline, projected if a
#' projection is current. The projection is determined in this order:
#' 1. An explicit `proj.out` argument
#' 2. The projection stored by the last [plot.indicatrix_list()] call
#'
#' `tissot_abline()` draws vertical and horizontal reference lines at a
#' given longitude/latitude in projected coordinates.
#'
#' @param proj.out target CRS. If `NULL`, uses the last plot projection
#'   (from `plot.indicatrix_list()`) or draws in lon/lat.
#' @param add logical; add to existing plot (default `TRUE`) or create new
#' @param ... graphical parameters passed to [graphics::lines()] (if adding)
#'   or [graphics::plot()] (if creating new)
#' @return `tissot_map()` invisibly returns the (projected) world coastline matrix
#' @seealso [plot.indicatrix_list()], [tissot()]
#' @export
#' @examples
#' r <- tissot(cbind(seq(-150, 150, by = 30), 0), proj.out = "+proj=robin")
#' ii <- indicatrix(r)
#' plot(ii, scale = 6e5, add = FALSE)
#' tissot_map()
tissot_map <- function(..., proj.out = NULL, add = TRUE) {
  props <- list(...)

  if (is.null(proj.out)) {
    proj.out <- getOption("tissot.last.plot.proj")
  }

  w <- if (is.null(proj.out)) {
    world
  } else {
    gdalraster::transform_xy(world, "OGC:CRS84", proj.out)
  }

  if (is.null(props$col)) props$col <- grDevices::rgb(.7, .7, .7)
  if (add) {
    props$x <- w
    do.call(graphics::lines, props)
  } else {
    props$x <- w
    if (is.null(props$pch)) props$pch <- "."
    do.call(graphics::plot, props)
    if (!is.null(proj.out)) {
      options(tissot.last.plot.proj = proj.out)
    }
  }
  invisible(w)
}


#' @rdname tissot_map
#' @param lambda longitude at which to draw a vertical line
#' @param phi latitude at which to draw a horizontal line
#' @param proj.in source CRS for lambda/phi (default "EPSG:4326")
#' @return `tissot_abline()` is called for its side effect
#' @export
#' @importFrom graphics abline
tissot_abline <- function(lambda, phi = NULL, ..., proj.in = NULL,
                          proj.out = NULL) {
  xy <- do.call(cbind, grDevices::xy.coords(lambda, phi)[1:2])

  if (is.null(proj.out)) {
    proj.out <- getOption("tissot.last.plot.proj")
  }
  if (!is.null(proj.out)) {
    if (is.null(proj.in)) proj.in <- "EPSG:4326"
    xy <- gdalraster::transform_xy(xy, proj.in, proj.out)
  }
  graphics::abline(v = xy[, 1L], h = xy[, 2L], ...)
}


#' Get or set the current plot projection
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' These functions access a global option. Prefer passing `proj.out`
#' explicitly to [tissot_map()] and [tissot_abline()], or let
#' [plot.indicatrix_list()] set the projection automatically.
#'
#' @return `tissot_get_proj()` returns the current projection string or `NULL`
#' @export
tissot_get_proj <- function() {
  getOption("tissot.last.plot.proj")
}

#' @rdname tissot_get_proj
#' @param proj.out projection CRS string
#' @export
tissot_set_proj <- function(proj.out) {
  options(tissot.last.plot.proj = proj.out)
}
