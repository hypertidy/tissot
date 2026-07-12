# Compute distortion surfaces on a regular grid

Create a raster-like grid of Tissot distortion metrics in the target
projection. The grid is laid out in projected coordinates,
back-projected to longitude/latitude, and passed through
[`tissot()`](https://hypertidy.github.io/tissot/reference/tissot.md).
Returns a list of matrices compatible with
[`image()`](https://rdrr.io/r/graphics/image.html) and easily converted
to raster objects (e.g. via `terra::rast()`).

## Usage

``` r
tissot_raster(
  target,
  extent = NULL,
  radius = 2e+07,
  nx = 100L,
  ny = NULL,
  metrics = c("scale_area", "angle_deformation", "scale_h", "scale_k"),
  ...,
  source = "EPSG:4326"
)

# S3 method for class 'tissot_raster'
image(x, metric = NULL, col = NULL, asp = 1, main = NULL, ...)

# S3 method for class 'tissot_raster'
plot(x, ...)
```

## Arguments

- target:

  target projection CRS string

- extent:

  numeric length 4: `c(xmin, xmax, ymin, ymax)` in projected
  coordinates. If `NULL` (default), estimated from the global lonlat
  bounding box (clamped by `radius`).

- radius:

  maximum half-width of the grid extent in projected units (default
  `2e7`). The auto-detected extent is clamped so that no dimension
  exceeds `2 * radius`. Prevents blow-up for projections like
  stereographic or gnomonic that map the whole globe to huge
  coordinates. Ignored if `extent` is supplied.

- nx:

  integer; number of grid cells in x (default 100)

- ny:

  integer; number of grid cells in y. If `NULL`, set to maintain
  approximately square cells.

- metrics:

  character vector of column names from
  [`tissot()`](https://hypertidy.github.io/tissot/reference/tissot.md)
  output to include as layers. Default includes areal scale, angular
  deformation, and meridional/parallel scale factors.

- ...:

  passed to [`image()`](https://rdrr.io/r/graphics/image.html)

- source:

  source CRS (default `"EPSG:4326"`)

- x:

  a `tissot_raster` (from `tissot_raster()`)

- metric:

  which metric to plot (default: first available)

- col:

  colour palette (default: `hcl.colors(64, "YlOrRd", rev = TRUE)`)

- asp:

  aspect ratio (default 1)

- main:

  plot title (default: metric name)

## Value

A `tissot_raster` list with components:

- `x`:

  numeric vector of x (easting) cell centers

- `y`:

  numeric vector of y (northing) cell centers

- `z`:

  named list of `ny x nx` matrices, one per metric

- `target`:

  the target CRS string

- `extent`:

  the `c(xmin, xmax, ymin, ymax)` used

The object has an [`image()`](https://rdrr.io/r/graphics/image.html)
method for quick plotting.

## Details

If `extent` is `NULL`, the function samples a dense set of boundary
coordinates spanning \\\[-180, 180\] \times \[-85, 85\]\\
(longitude/latitude), projects them into the target CRS, and uses the
resulting bounding range (padded by 5%) as the grid extent.

## See also

[`tissot()`](https://hypertidy.github.io/tissot/reference/tissot.md),
`plot.tissot_raster()`, `image.tissot_raster()`

## Examples

``` r
tr <- tissot_raster("+proj=robin", nx = 60)
image(tr)


## specific metric
image(tr, metric = "angle_deformation")


## with coastline overlay
image(tr)
tissot_map()


## polar stereo — use radius to cap extent
tp <- tissot_raster("+proj=stere +lat_0=-90", radius = 5e6, nx = 60)
image(tp)
tissot_map()
```
