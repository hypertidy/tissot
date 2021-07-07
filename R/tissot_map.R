#' Get last plot projection
#'
#'
#' 'tissot_map()' will add the [world] coastline to any map.
#'
#' 'tissot_get_proj()' When the indicatrix is plotted it registers its projection. This string
#' can be obtained with this getter function.
#'
#' @param add logical, default 'TRUE' add to existing plot or create new
#' @param ... graphical parameters for [lines()] if 'add = TRUE', or for [plot()] if 'add = FALSE'
#' @return 'tissot_map()' returns the internal world map data (projected if one is current) as a matrix
#' @export
tissot_map <- function(..., add = TRUE) {
  ## we might catch list(...) and see if something can project them? (one day)
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
#' @name tissot_abline
#' @export
#' @return 'tissot_abline()' called for its side effect of drawing on the plot
tissot_abline <- function(lambda, phi  = NULL, ..., proj.in = NULL) {
  xy <- do.call(cbind, xy.coords(lambda, phi)[1:2])

  target <- tissot_get_proj()
  if (!is.null(target)) {
    xy <- .prj(xy, target, proj.in = proj.in)
  }
  abline(v = xy[,1L], h = xy[,2L])
}
#' @name tissot_map
#' @return 'tissot_get_proj()' returns the value of the current projection, or NULL
#' @export
tissot_get_proj <- function() {
  getOption("tissot.last.plot.proj")
}

.project_world <- function() {
  target <- tissot_get_proj()
  if (is.null(target)) return(world)  ## no projection in effect
  .prj(world, target, proj.in = "OGC:CRS84")
}
