# tissot 0.2.0

* Switched projection engine from `reproj` to `gdalraster`. All coordinate
  transformations now use `gdalraster::transform_xy()`, supporting any CRS
  that GDAL/PROJ can handle.

* `tissot()` computes the Jacobian via a single batched projection call
  (3N points in one shot) with fully vectorized distortion metrics — no

  per-point loop.

* `tissot()` now returns a `tissot_tbl` subclass carrying `proj.in` and
  `proj.out` as attributes. New `print()` and `summary()` methods give a
  quick overview of distortion statistics.

* `indicatrix()` accepts either a `tissot_tbl` (projection extracted
  automatically from attributes) or raw lon/lat with explicit `proj.out`.
  The `scale` and `n` parameters moved to the plot methods where they belong.

* `indicatrix()` now returns an `indicatrix_list` class with `[`, `print`,
  and `length` methods, replacing the old unclassed list of `indicatrixes`.

* Rich plotting for `plot.indicatrix()`: reference unit circle
  (`show.circle`), lambda direction line in red and phi direction line in
  blue (`show.axes`), drawn behind and over the filled ellipse.

* `plot.indicatrix_list()` gains `fill.by` for colour-coded distortion
  (e.g. `fill.by = "scale.area"`) using `hcl.colors()` with no extra

  dependencies. Reference circles are shown by default to make distortion
  visible at map scale.

* `tissot_map()` and `tissot_abline()` gain an explicit `proj.out`
  argument. `tissot_get_proj()` and `tissot_set_proj()` are retained but
  soft-deprecated in favour of the projection carried by `indicatrix_list`
  attributes.

* Added comprehensive test suite via `testthat`: conformal, equal-area,
  and equidistant property checks, vectorization consistency, matrix input,
  class structure, polar edge cases, and method output.

* README rewritten to showcase the new API with projection comparisons,
  colour-coded distortion, and rich single-indicatrix plots.

# tissot 0.1.0

* New functions `tissot_map()` and `tissot_abline()` to easily add contextual
  world map data and location lines.
* Improved ability to build more than one tissot and to easily plot multiple
  indicatrixes.
* New function `tissot_get_proj()` with no arguments, just gets the last value
  of a plot projection registered by plot.indicatrixes.

# tissot 0.0.1

* Removed some old stuff.
* Using PROJ and libproj, dev only versions.
