# Unproject coordinates to geographic (lon/lat) CRS

A convenience wrapper around
[`PROJ::proj_trans()`](https://hypertidy.github.io/PROJ/reference/proj_trans.html)
for converting projected coordinates to geographic. Useful for
generating regular grids in a projected CRS to feed to
[`tissot()`](https://hypertidy.github.io/tissot/reference/tissot.md),
which requires lon/lat input.

## Usage

``` r
tissot_unproject(x, target = "EPSG:4326", ..., source = NULL)
```

## Arguments

- x:

  input coordinates — any xy-ish object: a two-column matrix,
  data.frame, tibble, list with `x`/`y` or `lon`/`lat` components, or a
  length-2 numeric vector for a single point

- target:

  target CRS string (default `"EPSG:4326"`). Must be geographic.

- ...:

  ignored

- source:

  source CRS string (required). Must be a projected CRS.

## Value

A two-column matrix of longitude and latitude values.

## See also

[`tissot()`](https://hypertidy.github.io/tissot/reference/tissot.md),
[`PROJ::proj_trans()`](https://hypertidy.github.io/PROJ/reference/proj_trans.html)

## Examples

``` r
## regular grid in UTM zone 55S, unprojected to lon/lat for tissot()
xy <- expand.grid(x = seq(4e5, 6e5, length.out = 5),
                  y = seq(5200000, 5400000, length.out = 4))
ll <- tissot_unproject(xy, source = "EPSG:32755")
tissot(ll, "+proj=utm +zone=55 +south")
#> Tissot indicatrix: 20 points, +proj=utm +zone=55 +south
#> # A tibble: 20 × 14
#>        x     y dx_dlam   dy_dlam  dx_dphi dy_dphi scale_h scale_k scale_omega
#>    <dbl> <dbl>   <dbl>     <dbl>    <dbl>   <dbl>   <dbl>   <dbl>       <dbl>
#>  1  146. -43.3   0.728  1.08e- 2 -0.0147    0.998   1.000   1.000 0          
#>  2  146. -43.4   0.728  5.38e- 3 -0.00737   0.998   1.000   1.000 0          
#>  3  147  -43.4   0.728  0         0         0.998   1.000   1.000 0          
#>  4  148. -43.4   0.728 -5.38e- 3  0.00737   0.998   1.000   1.000 0          
#>  5  148. -43.3   0.728 -1.08e- 2  0.0147    0.998   1.000   1.000 0          
#>  6  146. -42.7   0.735  1.06e- 2 -0.0144    0.998   1.000   1.000 0          
#>  7  146. -42.8   0.735  5.32e- 3 -0.00722   0.998   1.000   1.000 0.000000854
#>  8  147  -42.8   0.735  2.78e-12  0         0.998   1.000   1.000 0          
#>  9  148. -42.8   0.735 -5.32e- 3  0.00722   0.998   1.000   1.000 0.00000148 
#> 10  148. -42.7   0.735 -1.06e- 2  0.0144    0.998   1.000   1.000 0          
#> 11  146. -42.1   0.742  1.05e- 2 -0.0141    0.997   1.000   1.000 0.000000854
#> 12  146. -42.2   0.742  5.26e- 3 -0.00707   0.997   1.000   1.000 0          
#> 13  147  -42.2   0.742  0         0         0.997   1.000   1.000 0.00000121 
#> 14  148. -42.2   0.742 -5.26e- 3  0.00707   0.997   1.000   1.000 0          
#> 15  148. -42.1   0.742 -1.05e- 2  0.0141    0.997   1.000   1.000 0          
#> 16  146. -41.5   0.749  1.04e- 2 -0.0138    0.997   1.000   1.000 0          
#> 17  146. -41.6   0.749  5.20e- 3 -0.00692   0.997   1.000   1.000 0.00000121 
#> 18  147  -41.6   0.749  2.78e-12  0         0.997   1.000   1.000 0          
#> 19  148. -41.6   0.749 -5.20e- 3  0.00692   0.997   1.000   1.000 0.00000148 
#> 20  148. -41.5   0.749 -1.04e- 2  0.0138    0.997   1.000   1.000 0.000000854
#> # ℹ 5 more variables: scale_a <dbl>, scale_b <dbl>, scale_area <dbl>,
#> #   angle_deformation <dbl>, convergence <dbl>
```
