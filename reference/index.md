# Package index

## Distortion metrics

Compute Tissot indicatrix properties at arbitrary points or over a
regular projected grid. Returns scale factors, angular deformation,
areal distortion, and convergence.

- [`tissot()`](https://hypertidy.github.io/tissot/reference/tissot.md) :
  Compute Tissot indicatrix properties
- [`tissot_raster()`](https://hypertidy.github.io/tissot/reference/tissot_raster.md)
  [`image(`*`<tissot_raster>`*`)`](https://hypertidy.github.io/tissot/reference/tissot_raster.md)
  [`plot(`*`<tissot_raster>`*`)`](https://hypertidy.github.io/tissot/reference/tissot_raster.md)
  : Compute distortion surfaces on a regular grid

## Indicatrix objects and plotting

Build indicatrix objects from distortion results, extract raw ellipse
coordinates, and draw single or multiple indicatrixes on a map.

- [`indicatrix()`](https://hypertidy.github.io/tissot/reference/indicatrix.md)
  : Build Tissot indicatrix objects
- [`` `[`( ``*`<indicatrix_list>`*`)`](https://hypertidy.github.io/tissot/reference/sub-.indicatrix_list.md)
  : Subset an indicatrix_list
- [`ti_ellipse()`](https://hypertidy.github.io/tissot/reference/ti_ellipse.md)
  : Generate ellipse coordinates for an indicatrix
- [`plot(`*`<indicatrix>`*`)`](https://hypertidy.github.io/tissot/reference/plot.indicatrix.md)
  : Plot an indicatrix
- [`plot(`*`<indicatrix_list>`*`)`](https://hypertidy.github.io/tissot/reference/plot.indicatrix_list.md)
  : Plot a list of indicatrixes

## Raster distortion plots

Visualize a `tissot_raster` distortion surface with
[`image()`](https://rdrr.io/r/graphics/image.html) or
[`plot()`](https://rdrr.io/r/graphics/plot.default.html).

- [`tissot_raster()`](https://hypertidy.github.io/tissot/reference/tissot_raster.md)
  [`image(`*`<tissot_raster>`*`)`](https://hypertidy.github.io/tissot/reference/tissot_raster.md)
  [`plot(`*`<tissot_raster>`*`)`](https://hypertidy.github.io/tissot/reference/tissot_raster.md)
  : Compute distortion surfaces on a regular grid

## Map utilities

Add a world coastline overlay, reference ablines at given lon/lat
positions, or convert projected coordinates back to geographic for use
with
[`tissot()`](https://hypertidy.github.io/tissot/reference/tissot.md).

- [`tissot_map()`](https://hypertidy.github.io/tissot/reference/tissot_map.md)
  [`tissot_abline()`](https://hypertidy.github.io/tissot/reference/tissot_map.md)
  : Plot world coastline on a projected map
- [`tissot_unproject()`](https://hypertidy.github.io/tissot/reference/tissot_unproject.md)
  : Unproject coordinates to geographic (lon/lat) CRS

## Data

Bundled datasets used in examples.

- [`world`](https://hypertidy.github.io/tissot/reference/world.md) :
  World coastline

## Deprecated

Superseded helpers for accessing the package-level projection state.
Prefer passing `target` explicitly to
[`tissot_map()`](https://hypertidy.github.io/tissot/reference/tissot_map.md)
and
[`tissot_abline()`](https://hypertidy.github.io/tissot/reference/tissot_map.md).

- [`tissot_get_proj()`](https://hypertidy.github.io/tissot/reference/tissot_get_proj.md)
  [`tissot_set_proj()`](https://hypertidy.github.io/tissot/reference/tissot_get_proj.md)
  : Get or set the current plot projection
