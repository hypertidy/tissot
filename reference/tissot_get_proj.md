# Get or set the current plot projection

**Deprecated.** Prefer passing `target` explicitly.

These functions access the package-level projection state. Prefer
passing `target` explicitly to
[`tissot_map()`](https://hypertidy.github.io/tissot/reference/tissot_map.md)
and
[`tissot_abline()`](https://hypertidy.github.io/tissot/reference/tissot_map.md),
or let
[`plot.indicatrix_list()`](https://hypertidy.github.io/tissot/reference/plot.indicatrix_list.md)
or
[`image.tissot_raster()`](https://hypertidy.github.io/tissot/reference/tissot_raster.md)
set it automatically.

## Usage

``` r
tissot_get_proj()

tissot_set_proj(target)
```

## Arguments

- target:

  projection CRS string

## Value

`tissot_get_proj()` returns the current projection string or `NULL`

## Examples

``` r
tissot_set_proj("+proj=robin")
tissot_get_proj()
#> [1] "+proj=robin"
tissot_set_proj(NULL)  # reset
```
