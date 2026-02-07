#' @keywords internal
#' @aliases tissot-package
#' @importFrom gdalraster transform_xy
#' @importFrom graphics lines polygon abline
#' @importFrom grDevices adjustcolor hcl.colors rgb
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL

globalVariables("world")

#' World coastline
#'
#' A matrix of longitude/latitude coordinates representing a simplified
#' world coastline, derived from the maps package. Continent boundaries
#' are separated by `NA` rows. Longitudes are constrained to
#' `abs(lon) <= 180`.
#'
#' @format A two-column numeric matrix (longitude, latitude)
#' @docType data
#' @name world
"world"
