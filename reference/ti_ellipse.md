# Generate ellipse coordinates for an indicatrix

Generate ellipse coordinates for an indicatrix

## Usage

``` r
ti_ellipse(x, scale = 1e+05, n = 72, ...)
```

## Arguments

- x:

  an indicatrix object

- scale:

  scaling factor for the ellipse size in projected units

- n:

  number of points on the ellipse

- ...:

  ignored

## Value

A two-column matrix of projected coordinates tracing the ellipse

## See also

[`indicatrix()`](https://hypertidy.github.io/tissot/reference/indicatrix.md),
[`plot.indicatrix()`](https://hypertidy.github.io/tissot/reference/plot.indicatrix.md)

## Examples

``` r
ii <- indicatrix(c(0, 45), "+proj=robin")
ell <- ti_ellipse(ii[[1]], scale = 5e5)
head(ell)
#>          [,1]    [,2]
#> [1,] 380302.5 4805074
#> [2,] 378855.3 4845730
#> [3,] 374524.8 4886078
#> [4,] 367344.0 4925808
#> [5,] 357367.4 4964620
#> [6,] 344671.1 5002218

## draw on a map
tissot_map(target = "+proj=robin", add = FALSE)
lines(ell)
```
