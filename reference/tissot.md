# Compute Tissot indicatrix properties

Compute the Tissot indicatrix at given longitude/latitude locations for
a map projection. Returns scale factors, angular deformation,
convergence, and related distortion properties.

## Usage

``` r
tissot(
  x,
  target,
  ...,
  source = "EPSG:4326",
  method = c("proj", "finitediff"),
  A = 6378137,
  f.inv = 298.257223563,
  dx = 1e-04
)
```

## Arguments

- x:

  input coordinates — any xy-ish object: a two-column matrix,
  data.frame, tibble, list with `x`/`y` or `lon`/`lat` components, or a
  length-2 numeric vector for a single point

- target:

  target projection CRS string (required)

- ...:

  ignored

- source:

  source CRS (default `"EPSG:4326"`)

- method:

  computation method: `"proj"` (default) uses
  [`PROJ::proj_factors()`](https://hypertidy.github.io/PROJ/reference/proj_factors.html);
  `"finitediff"` uses a finite-difference Jacobian (Snyder 1987)

- A:

  semi-major axis of the ellipsoid (default WGS84;
  `method = "finitediff"` only)

- f.inv:

  inverse flattening (default WGS84; `method = "finitediff"` only)

- dx:

  finite difference step in degrees (default 1e-4;
  `method = "finitediff"` only)

## Value

A `tissot_tbl` tibble with columns: `x` (lon), `y` (lat), `dx_dlam`,
`dy_dlam`, `dx_dphi`, `dy_dphi`, `scale_h`, `scale_k`, `scale_omega`,
`scale_a`, `scale_b`, `scale_area`, `angle_deformation`, `convergence`.
The `source` and `target` CRS strings are stored as attributes.

## Details

By default (`method = "proj"`) distortion measures are computed via
[`PROJ::proj_factors()`](https://hypertidy.github.io/PROJ/reference/proj_factors.html),
which calls the PROJ C library directly and is more accurate than finite
differences (particularly for tabular or piecewise-defined projections).
The original finite-difference path based on Snyder (1987) — inspired by
Bill Huber's formulation at <https://gis.stackexchange.com/a/5075/482> —
is preserved as `method = "finitediff"`.

`x` is assumed to contain longitude/latitude values; the default source
CRS is `"EPSG:4326"`. Set `source` for a different geographic CRS.

## See also

[`indicatrix()`](https://hypertidy.github.io/tissot/reference/indicatrix.md),
[`tissot_map()`](https://hypertidy.github.io/tissot/reference/tissot_map.md),
[`tissot_unproject()`](https://hypertidy.github.io/tissot/reference/tissot_unproject.md),
[`PROJ::proj_factors()`](https://hypertidy.github.io/PROJ/reference/proj_factors.html),
[`PROJ::proj_trans()`](https://hypertidy.github.io/PROJ/reference/proj_trans.html)

## Examples

``` r
tissot(c(0, 45), "+proj=robin")
#> Tissot indicatrix: 1 point, +proj=robin
#> # A tibble: 1 × 14
#>       x     y dx_dlam  dy_dlam dx_dphi dy_dphi scale_h scale_k scale_omega
#>   <dbl> <dbl>   <dbl>    <dbl>   <dbl>   <dbl>   <dbl>   <dbl>       <dbl>
#> 1     0    45   0.761 2.78e-12       0   0.933   0.933    1.08        8.15
#> # ℹ 5 more variables: scale_a <dbl>, scale_b <dbl>, scale_area <dbl>,
#> #   angle_deformation <dbl>, convergence <dbl>
tissot(cbind(seq(-180, 180, by = 30), 0), "+proj=robin")
#> Tissot indicatrix: 13 points, +proj=robin
#> # A tibble: 13 × 14
#>        x     y dx_dlam dy_dlam dx_dphi dy_dphi scale_h scale_k scale_omega
#>    <dbl> <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>       <dbl>
#>  1  -180     0   0.849       0       0   0.961   0.961   0.849        7.10
#>  2  -150     0   0.849       0       0   0.961   0.961   0.849        7.10
#>  3  -120     0   0.849       0       0   0.961   0.961   0.849        7.10
#>  4   -90     0   0.849       0       0   0.961   0.961   0.849        7.10
#>  5   -60     0   0.849       0       0   0.961   0.961   0.849        7.10
#>  6   -30     0   0.849       0       0   0.961   0.961   0.849        7.10
#>  7     0     0   0.849       0       0   0.961   0.961   0.849        7.10
#>  8    30     0   0.849       0       0   0.961   0.961   0.849        7.10
#>  9    60     0   0.849       0       0   0.961   0.961   0.849        7.10
#> 10    90     0   0.849       0       0   0.961   0.961   0.849        7.10
#> 11   120     0   0.849       0       0   0.961   0.961   0.849        7.10
#> 12   150     0   0.849       0       0   0.961   0.961   0.849        7.10
#> 13   180     0   0.849       0       0   0.961   0.961   0.849        7.10
#> # ℹ 5 more variables: scale_a <dbl>, scale_b <dbl>, scale_area <dbl>,
#> #   angle_deformation <dbl>, convergence <dbl>

## compare methods
tissot(c(0, 45), "+proj=robin", method = "proj")
#> Tissot indicatrix: 1 point, +proj=robin
#> # A tibble: 1 × 14
#>       x     y dx_dlam  dy_dlam dx_dphi dy_dphi scale_h scale_k scale_omega
#>   <dbl> <dbl>   <dbl>    <dbl>   <dbl>   <dbl>   <dbl>   <dbl>       <dbl>
#> 1     0    45   0.761 2.78e-12       0   0.933   0.933    1.08        8.15
#> # ℹ 5 more variables: scale_a <dbl>, scale_b <dbl>, scale_area <dbl>,
#> #   angle_deformation <dbl>, convergence <dbl>
tissot(c(0, 45), "+proj=robin", method = "finitediff")
#> Tissot indicatrix: 1 point, +proj=robin
#> # A tibble: 1 × 14
#>       x     y dx_dlam dy_dlam dx_dphi dy_dphi scale_h scale_k scale_omega
#>   <dbl> <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>       <dbl>
#> 1     0    45    1.07       0       0   0.946   0.946    1.07        7.25
#> # ℹ 5 more variables: scale_a <dbl>, scale_b <dbl>, scale_area <dbl>,
#> #   angle_deformation <dbl>, convergence <dbl>
```
