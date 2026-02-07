# Tissot Package Refactor Plan — February 2026

**Branch:** `refactor-2026`

## ✅ Completed

### 1. Switched projection engine from `reproj` to `gdalraster`
- `DESCRIPTION`: `reproj` → `gdalraster` in Imports
- `R/tissot.R`: Complete rewrite of `tissot()` to batch all finite-difference
  points (base + λ-offset + φ-offset = 3N points) into a **single
  `gdalraster::transform_xy()` call**, then compute Jacobian and all Tissot
  properties with vectorized R (no loop).
- `R/tissot_map.R`: `reproj::reproj_xy` → `gdalraster::transform_xy`
- `R/tissot-package.R`: Updated import
- `NAMESPACE` + `man/`: Regenerated via `devtools::document()`
- Removed stale `tissot/DESCRIPTION` subdirectory copy

### 2. Smoke-tested the new engine
- Equal-area (LAEA): `scale.area = 1` ✓
- Conformal (Stereographic): `angle_deformation ≈ 0` (5e-05 numerical noise) ✓
- Robinson vectorized: 13 points computed in one shot ✓
- Indicatrix object structure correct ✓

### 3. ✅ Reconciled `indicatrix()` API with old usage patterns

**What was done:**

- Added `"tissot_tbl"` subclass to the tibble returned by `tissot()`, with
  `proj.in` and `proj.out` stored as attributes (not repeated per row)
- Added `print.tissot_tbl()` and `summary.tissot_tbl()` methods
- `indicatrix()` now accepts either:
  - A `tissot_tbl` object (detects via class; extracts proj from attributes)
  - Raw lon/lat + `proj.out` (direct path)
- `indicatrix()` returns an `indicatrix_list` (classed list) with `proj.in`
  and `proj.out` as attributes
- Added `[.indicatrix_list`, `print.indicatrix_list`, `length.indicatrix_list`
  methods for ergonomic subsetting
- `scale` and `n` are now parameters of `plot()`, not `indicatrix()` —
  these are display parameters, not geometric properties

### 4. ✅ Restored rich plotting

**What was done:**

- `plot.indicatrix()` now draws:
  - Reference unit circle (dashed white, behind ellipse) via `show.circle`
  - Filled indicatrix ellipse
  - Lambda direction line (red) — parallel scale, via `show.axes`
  - Phi direction line (blue) — meridional scale, via `show.axes`
- `plot.indicatrix_list()`:
  - Iterates all indicatrices with `add = TRUE`
  - Creates new plot from bounding box when `add = FALSE`
  - Auto-registers projection for `tissot_map()` integration
  - Colour-coded fill via `fill.by = "scale.area"` etc., using
    `hcl.colors()` (no extra deps)
- Internal helpers: `ti_circle()`, `ti_axes()`

### 5. ✅ Removed global projection state (mostly)

**What was done:**

- `tissot_map()` now accepts explicit `proj.out` argument
- `tissot_abline()` now accepts explicit `proj.out` argument
- Both fall back to `getOption("tissot.last.plot.proj")` if `proj.out = NULL`
  (set automatically by `plot.indicatrix_list()`)
- `tissot_get_proj()` and `tissot_set_proj()` retained but documented as
  deprecated — the `indicatrix_list` carries `proj.out` as an attribute,
  making explicit state the primary path

### 6. ✅ Set up tests (`testthat`)

**File:** `tests/testthat/test-tissot.R`

**Test cases implemented:**
- Conformal projection: `angle_deformation ≈ 0`, `scale.a ≈ scale.b`
- Equal-area projection: `scale.area ≈ 1`
- Plate Carrée at equator: `scale.h ≈ 1`, `scale.k ≈ 1`
- Vectorized consistency: multi-point matches single-point calls
- Matrix input: `tissot(cbind(lon, lat), ...)` same as `tissot(lon, lat, ...)`
- `tissot_tbl` class and attributes
- All expected columns present
- `indicatrix()` from `tissot_tbl` (item 3 verification)
- `indicatrix()` from raw lon/lat
- `indicatrix()` errors without `proj.out`
- `ti_ellipse()` returns closed polygon (first == last row)
- `indicatrix_list` subsetting preserves class and attributes
- Multiple conformal projections (Stereographic, Mercator, TM)
- Multiple equal-area projections (LAEA, Mollweide, Sinusoidal)
- Polar points don't produce NaN
- `scale.omega` equals `angle_deformation`
- `print` and `summary` methods produce expected output

### 7. ✅ Updated README.Rmd

**What was done:**
- Complete rewrite — removed all references to:
  - `indicatrix0()` (removed function)
  - `tissot:::.prj()` (removed internal)
  - `earthcircle::earthcircle()` (not a dependency)
  - `maps::map()` (not a dependency)
  - Old column names (`axes_*`, `lambda_d*`, `phi_d*`)
- New sections:
  - Quick start: `tissot()` → tibble of distortion properties
  - Column reference table
  - Plotting indicatrices with coastline overlay
  - Projection comparison (Robinson / Mollweide / Mercator side by side)
  - Colour-coded distortion (`fill.by`)
  - Rich single-indicatrix plots (axes + reference circle)
  - Polar projections
  - Distortion summary

### 8. ✅ Clean up man pages

- Stale `.Rd` files for `indicatrix0`, `indicatrixes` already cleaned
- NAMESPACE updated with all new S3 methods
- `@seealso` cross-references between functions
- `@param` docs match current signatures
- **Note:** run `devtools::document()` to regenerate `man/` from updated roxygen

---

## 9. Future opportunities (not blocking this refactor)

- `tissot_grid(proj.out, extent, spacing)` — convenience for regular grids
- `tissot_graticule()` — draw projected graticule lines
- Natural Earth coastline to replace the `{maps}`-derived `world` data
- `geographiclib` (via `geoarea`) in Suggests for geodesic-exact validation
- Colour-coded distortion surface via `image()` / contour
- Projection comparison: `tissot_compare(proj1, proj2, ...)`
- Vignette: `vignette("tissot")` with mathematical background

## File inventory (current state on branch)

```
R/tissot.R          — Core engine + tissot_tbl class + indicatrix + indicatrix_list + plotting
R/tissot_map.R      — Map overlay, abline, deprecated get/set_proj
R/tissot-package.R  — Package docs, imports, world data docs
R/utils.R           — Internal %||% operator
DESCRIPTION         — v0.1.0, gdalraster, testthat in Suggests
NAMESPACE           — All S3 methods registered
tests/testthat.R    — Test runner
tests/testthat/test-tissot.R — Comprehensive test suite
README.Rmd          — Full rewrite for new API
data/world.rda      — Bundled coastline matrix (unchanged)
data-raw/DATASET.R  — (unchanged)
man/                — Needs regeneration via devtools::document()
```
