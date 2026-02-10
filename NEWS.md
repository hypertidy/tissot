# tissot 0.2.0

Major refactor — modernized engine, new API, and rich plotting.

## Breaking changes

* Projection engine switched from `reproj` to `gdalraster::transform_xy()`.
  All projection work is now a single batched GDAL call.

* API renamed: `proj.in`/`proj.out` → `source`/`target`, following reproj
  conventions. `target` is the second positional argument (required);
  `source` defaults to `"EPSG:4326"`.

* `tissot()` input is now any "xy-ish" object (matrix, data.frame, list,
  length-2 vector) via the internal `as_xy()` helper. The old
  `lambda`/`phi` positional arguments are replaced by `x`.

* `tissot()` returns a `tissot_tbl` (subclassed tibble) with `source` and
  `target` stored as attributes. Column names changed to: `x`, `y`,
  `dx_dlam`, `dy_dlam`, `dx_dphi`, `dy_dphi`, `scale.h`, `scale.k`,
  `scale.omega`, `scale.a`, `scale.b`, `scale.area`, `angle_deformation`,
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
    - Colour-coded fill via `fill.by` (e.g. `"scale.area"`,
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
