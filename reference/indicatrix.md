# Build Tissot indicatrix objects

Computes an `indicatrix_list` containing per-point distortion objects
suitable for plotting with
[`plot.indicatrix_list()`](https://hypertidy.github.io/tissot/reference/plot.indicatrix_list.md)
or querying directly.

## Usage

``` r
indicatrix(x, target = NULL, ..., source = "EPSG:4326")
```

## Arguments

- x:

  a `tissot_tbl`, or any xy-ish input (see
  [`tissot()`](https://hypertidy.github.io/tissot/reference/tissot.md))

- target:

  target projection CRS (extracted from `tissot_tbl` attributes if `x`
  is one; required otherwise)

- ...:

  passed to
  [`tissot()`](https://hypertidy.github.io/tissot/reference/tissot.md)

- source:

  source CRS (default `"EPSG:4326"`)

## Value

An `indicatrix_list` object (a list of `indicatrix` objects with
`source` and `target` stored as attributes)

## Details

`indicatrix()` accepts either:

- A `tissot_tbl` object (from
  [`tissot()`](https://hypertidy.github.io/tissot/reference/tissot.md))
  — projection is extracted from attributes, `target` is optional

- Any xy-ish input with an explicit `target`

## See also

[`tissot()`](https://hypertidy.github.io/tissot/reference/tissot.md),
[`plot.indicatrix()`](https://hypertidy.github.io/tissot/reference/plot.indicatrix.md),
[`plot.indicatrix_list()`](https://hypertidy.github.io/tissot/reference/plot.indicatrix_list.md),
[`ti_ellipse()`](https://hypertidy.github.io/tissot/reference/ti_ellipse.md)

## Examples

``` r
## From a tissot_tbl
r <- tissot(cbind(seq(-150, 150, by = 30), 0), "+proj=robin")
ii <- indicatrix(r)
plot(ii, scale = 5e5)


## From raw coordinates
ii2 <- indicatrix(c(0, 45), "+proj=stere +lat_0=90")
plot(ii2)
```
