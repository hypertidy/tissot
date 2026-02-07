# Null-coalescing operator (base R has this from 4.4.0 but
# we support R >= 3.6.0 per DESCRIPTION)
`%||%` <- function(x, y) if (is.null(x)) y else x


#' Resolve a logical-or-list graphical parameter
#'
#' If `x` is `TRUE`, return `defaults`. If `x` is a list, merge it over
#' `defaults` (user values win). If `x` is `FALSE`/`NULL`, return `NULL`.
#'
#' @param x `TRUE`, `FALSE`, or a named list of graphical parameters
#' @param defaults named list of default graphical parameters
#' @return A named list, or `NULL` if drawing is suppressed
#' @keywords internal
resolve_gpar <- function(x, defaults) {
  if (is.list(x)) {
    out <- defaults
    out[names(x)] <- x
    return(out)
  }
  if (isTRUE(x)) return(defaults)
  NULL
}
