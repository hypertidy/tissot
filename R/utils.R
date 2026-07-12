# Null-coalescing operator (base R has this from 4.4.0 but
# we support R >= 3.6.0 per DESCRIPTION)
`%||%` <- function(x, y) if (is.null(x)) y else x


## Package-level state: avoids options() side-effects that trigger linters.
## Access via .tissot_state$last_proj rather than options(tissot.last.plot.proj).
.tissot_state <- new.env(parent = emptyenv())
.tissot_state$last_proj <- NULL


## Internal CRS predicate: returns TRUE if the CRS is geographic (lon/lat).
## Uses PROJ::proj_crs_text() to get WKT2; geographic CRS WKT2 always
## begins with "GEOGCRS[" or "GEOGRAPHICCRS[".  Falls back to inspecting
## the PROJ string for +proj=longlat|latlong when WKT conversion returns NA
## (can happen for bare operation strings without +type=crs).
.srs_is_geographic <- function(crs) {
  wkt <- tryCatch(
    PROJ::proj_crs_text(crs, format = 0L),
    error = function(e) NA_character_
  )
  if (!is.na(wkt)) {
    grepl("^GEOG", wkt)
  } else {
    grepl("\\+proj=(longlat|latlong|lonlat)", crs)
  }
}


#' Coerce xy-ish input to a two-column matrix
#'
#' Accepts a two-column matrix, data.frame, tibble, list with x/y or
#' lon/lat components, or a length-2 numeric vector (single point).
#' Always returns an N x 2 numeric matrix.
#'
#' @param x an xy-ish object
#' @return A two-column numeric matrix
#' @keywords internal
#' @noRd
as_xy <- function(x) {
  if (is.matrix(x)) {
    stopifnot(ncol(x) >= 2L)
    return(matrix(as.double(x[, 1:2]), ncol = 2L))
  }
  if (is.data.frame(x)) {
    return(cbind(as.double(x[[1L]]), as.double(x[[2L]])))
  }
  if (is.list(x)) {
    nms <- names(x)
    if ("x" %in% nms && "y" %in% nms) {
      return(cbind(as.double(x$x), as.double(x$y)))
    }
    if ("lon" %in% nms && "lat" %in% nms) {
      return(cbind(as.double(x$lon), as.double(x$lat)))
    }
    if (length(x) >= 2L) {
      return(cbind(as.double(x[[1L]]), as.double(x[[2L]])))
    }
  }
  if (is.numeric(x) && length(x) == 2L) {
    return(matrix(as.double(x), nrow = 1L, ncol = 2L))
  }
  stop("Cannot coerce input to an xy matrix. Supply a two-column matrix, ",
       "data.frame, or list with x/y components.")
}


#' Resolve a logical-or-list graphical parameter
#'
#' If `x` is `TRUE`, return `defaults`. If `x` is a list, merge it over
#' `defaults` (user values win). If `x` is `FALSE`/`NULL`, return `NULL`.
#'
#' @param x `TRUE`, `FALSE`, or a named list of graphical parameters
#' @param defaults named list of default graphical parameters
#' @return A named list, or `NULL` if drawing is suppressed
#' @keywords internal
#' @noRd
resolve_gpar <- function(x, defaults) {
  if (is.list(x)) {
    out <- defaults
    out[names(x)] <- x
    return(out)
  }
  if (isTRUE(x)) return(defaults)
  NULL
}
