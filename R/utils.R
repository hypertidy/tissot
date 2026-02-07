# Null-coalescing operator (base R has this from 4.4.0 but
# we support R >= 2.10 per DESCRIPTION)
`%||%` <- function(x, y) if (is.null(x)) y else x
