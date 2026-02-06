#' Get last plot projection
#'
#'
#' 'tissot_map()' will add the [world] coastline to any map.
#'
#' 'tissot_get_proj()' When the indicatrix is plotted it registers its projection. This string
#' can be obtained with this getter function.
#'
#' 'tissot_abline()' will draw a vertical and horizontal line at a give longitude latitude (where
#' they intersect is the actual lon,lat location)
#'
#' @param add logical, default 'TRUE' add to existing plot or create new
#' @param ... graphical parameters for [lines()] if 'add = TRUE', or for [plot()] if 'add = FALSE'
#' @param lambda longitude at which to draw a vertical line
#' @param phi latitude at which to draw a horizontal line
#' @param proj.in projection for expert use
#' @return 'tissot_map()' returns the internal world map data (projected if one is current) as a matrix
#' @export
tissot_map <- function(..., add = TRUE) {
  props <- list(...)
  w <- .project_world()
  if (is.null(props$col)) props$col <- rgb(.7, .7, .7)
  if (add) {
    props$x <- w
    do.call(lines, props)
  } else {
    props$x <- w
    if (is.null(props$pch)) props$pch <- "."

    do.call(plot, props)
  }
  invisible(w)
}
#' @name tissot_map
#' @export
#' @importFrom graphics abline
#' @return 'tissot_abline()' called for its side effect of drawing on the plot
tissot_abline <- function(lambda, phi  = NULL, ..., proj.in = NULL) {
  xy <- do.call(cbind, xy.coords(lambda, phi)[1:2])

  target <- tissot_get_proj()
  if (!is.null(target)) {
    if (is.null(proj.in)) proj.in <- "EPSG:4326"
    xy <- gdalraster::transform_xy(xy, proj.in, target)
  }
  graphics::abline(v = xy[,1L], h = xy[,2L])
}
#' @name tissot_map
#' @return 'tissot_get_proj()' returns the value of the current projection, or NULL
#' @export
tissot_get_proj <- function() {
  getOption("tissot.last.plot.proj")
}

.project_world <- function() {
  target <- tissot_get_proj()
  if (is.null(target)) return(world)
  gdalraster::transform_xy(world, "OGC:CRS84", target)
}
