# tissot 0.3.0

## New features

* `tissot()` gains a `method` argument (default `"proj"`). The new `"proj"`
  method delegates to `PROJ::proj_factors()`, calling the PROJ C library
  directly for exact distortion factors rather than approximating them via
  finite differences. This is more accurate, especially for tabular or
  piecewise projections (e.g. Robinson). The original finite-difference
  computation (Snyder 1987), inspired by Bill Huber's formulation at
  <https://gis.stackexchange.com/a/5075/482>, is preserved as
  `method = "finitediff"`. The `PROJ` package is now the sole coordinate
  transformation dependency; `gdalraster` has been removed.

* `tissot_raster()` computes distortion surfaces on a regular projected grid,
  returning a `tissot_raster` object with `image()` and `plot()` methods.
  Extent is auto-detected from the global lon/lat bounding box (clamped by
  `radius`) or supplied explicitly. Supports any metric from [tissot()] output.

## Dependency change

* `gdalraster` removed from `Imports`. All coordinate transformation now goes
  through `PROJ::proj_trans()` and `PROJ::proj_crs_text()`. CRS validation
  (`srs_is_geographic`) is handled by an internal helper using
  `PROJ::proj_crs_text()`.

## Bug fixes

* Fixed default `metrics` in `tissot_raster()`: names now use underscores
  (`"scale_area"`, `"scale_h"`, `"scale_k"`) matching the actual `tissot_tbl`
  column names (previously used dots, causing silently all-NA output).

# tissot 0.2.0

Major refactor — modernized engine, new API, and rich plotting.

## Breaking changes

* Projection engine switched from `reproj` to `transform_xy()` from gdalraster.
  All projection work is now a single batched GDAL call.

* API renamed: `proj.in`/`proj.out` → `source`/`target`, following reproj
  conventions. `target` is the second positional argument (required);
  `source` defaults to `"EPSG:4326"`.

* `tissot()` input is now any "xy-ish" object (matrix, data.frame, list,
  length-2 vector) via the internal `as_xy()` helper. The old
  `lambda`/`phi` positional arguments are replaced by `x`.

* `tissot()` returns a `tissot_tbl` (subclassed tibble) with `source` and
  `target` stored as attributes. Column names changed to: `x`, `y`,
  `dx_dlam`, `dy_dlam`, `dx_dphi`, `dy_dphi`, `scale_h`, `scale_k`,
  `scale_omega`, `scale_a`, `scale_b`, `scale_area`, `angle_deformation`,
  `convergence`.

* `indicatrix()` returns an `indicatrix_list` (replaces `indicatrixes`
  class). Accepts `tissot_tbl` or raw coordinates + explicit `target`.

* Removed: `indicatrix0()`, `tissot0()`, `.prj()`, and the internal
  `numericDeriv` path. The Jacobian is now computed directly via
  finite differences.

## New features

* `plot.indicatrix()` and `plot.indicatrix_list()` with:
    - Reference unit circle overlay (`show.circle`)
    - Lambda/phi direction axes (`show.axes`)
    - `show.circle` and `show.axes` accept `TRUE`, `FALSE`, or a named list
      of graphical parameters for full customization
    - Colour-coded fill via `fill.by` (e.g. `"scale_area"`,
      `"angle_deformation"`)

* `print()` and `summary()` methods for `tissot_tbl`.

* `print()`, `length()`, and `[` methods for `indicatrix_list`.

* `tissot_map()` and `tissot_abline()` accept explicit `target` argument.

* `as_xy()` internal helper coerces diverse coordinate inputs.

* `resolve_gpar()` internal helper for the logical-or-list graphical
  parameter pattern.

## Internal

* Fully vectorized Jacobian computation — no per-point loop.
* Computation uses 3N batched `transform_xy()` call.
* Test suite added (testthat).
