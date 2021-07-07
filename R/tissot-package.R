#' @keywords internal
#' @aliases tissot-package
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL

globalVariables("world")
#' world coastline
#'
#' A modified matrix version of data from the maps package.
#'
#' Basically longitudes have been smooshed to 'abs(lon) < 180'
#' @docType data
#' @name world
'world'
