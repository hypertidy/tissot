# Plot an indicatrix

Draws a single Tissot indicatrix ellipse on the current plot. The
ellipse shows the distortion of a unit circle under the map projection.
Optional overlays include a reference unit circle, and lambda/phi
direction axes.

## Usage

``` r
# S3 method for class 'indicatrix'
plot(
  x,
  scale = 1e+05,
  n = 72,
  col = "#FF990055",
  border = "black",
  add = TRUE,
  show.axes = TRUE,
  show.circle = TRUE,
  ...
)
```

## Arguments

- x:

  an `indicatrix` object (from
  [`indicatrix()`](https://hypertidy.github.io/tissot/reference/indicatrix.md))

- scale:

  scaling factor for the ellipse size in projected units

- n:

  number of points on the ellipse

- col:

  fill colour for the ellipse

- border:

  border colour

- add:

  logical; add to existing plot (default `TRUE`)?

- show.axes:

  `TRUE`, `FALSE`, or a named list of graphical parameters for the
  direction lines. Defaults:
  `list(col.lambda = "red", col.phi = "blue", lwd = 1.5)`.

- show.circle:

  `TRUE`, `FALSE`, or a named list of graphical parameters for the
  reference circle. Defaults:
  `list(col = adjustcolor("white", alpha.f = 1), border = "grey70", lwd = 4.5, lty = 2)`.

- ...:

  passed to
  [`graphics::polygon()`](https://rdrr.io/r/graphics/polygon.html)

## Value

The input `x`, invisibly.

## Details

`show.circle` and `show.axes` accept `TRUE` (use defaults), `FALSE`
(hide), or a named list of graphical parameters to override defaults.
For example, `show.circle = list(border = "blue", lty = 3)`.

## See also

[`indicatrix()`](https://hypertidy.github.io/tissot/reference/indicatrix.md),
[`plot.indicatrix_list()`](https://hypertidy.github.io/tissot/reference/plot.indicatrix_list.md),
[`ti_ellipse()`](https://hypertidy.github.io/tissot/reference/ti_ellipse.md)
