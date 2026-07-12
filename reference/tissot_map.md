# Plot world coastline on a projected map

`tissot_map()` draws the bundled
[world](https://hypertidy.github.io/tissot/reference/world.md)
coastline, projected if a projection is current. The projection is
determined in this order:

1.  An explicit `target` argument

2.  The projection recorded by the most recent call to
    [`plot.indicatrix_list()`](https://hypertidy.github.io/tissot/reference/plot.indicatrix_list.md),
    [`image.tissot_raster()`](https://hypertidy.github.io/tissot/reference/tissot_raster.md),
    or
    [`tissot_raster()`](https://hypertidy.github.io/tissot/reference/tissot_raster.md)

## Usage

``` r
tissot_map(..., target = NULL, add = TRUE)

tissot_abline(x, y = NULL, ..., source = "EPSG:4326", target = NULL)
```

## Arguments

- ...:

  graphical parameters passed to
  [`graphics::lines()`](https://rdrr.io/r/graphics/lines.html) (if
  adding) or
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html) (if
  creating new)

- target:

  target CRS. If `NULL`, uses the last plot projection (from
  [`plot.indicatrix_list()`](https://hypertidy.github.io/tissot/reference/plot.indicatrix_list.md))
  or draws in lon/lat.

- add:

  logical; add to existing plot (default `TRUE`) or create new

- x:

  longitude values (or any xy-ish input; see
  [`tissot()`](https://hypertidy.github.io/tissot/reference/tissot.md))

- y:

  latitude values (ignored if `x` is a matrix)

- source:

  source CRS for the coordinates (default `"EPSG:4326"`)

## Value

`tissot_map()` invisibly returns the (projected) world coastline matrix

`tissot_abline()` is called for its side effect

## Details

`tissot_abline()` draws vertical and horizontal reference lines at a
given longitude/latitude in projected coordinates.

## See also

[`plot.indicatrix_list()`](https://hypertidy.github.io/tissot/reference/plot.indicatrix_list.md),
[`tissot()`](https://hypertidy.github.io/tissot/reference/tissot.md)

## Examples

``` r
r <- tissot(cbind(seq(-150, 150, by = 30), 0), "+proj=robin")
ii <- indicatrix(r)
plot(ii, scale = 6e5, add = FALSE)
tissot_map()
```
