# Tissot Package Refactor Plan — February 2026

## Branch: `refactor-2026`

## Completed

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

---

## Remaining Work

### 3. Reconcile `indicatrix()` API with old usage patterns

**Problem:** The old API was:
```r
r <- tissot(xy, proj.in = ..., proj.out = ...)
ii <- indicatrix(r, scale = 4e5, n = 36)   # took a tissot tibble
plot(ii, add = TRUE)                         # dispatched to plot.indicatrixes
```

The current rewrite changed `indicatrix()` to take lon/lat + `proj.out` directly
and returns a plain list (no class), so:

- `indicatrix(tibble)` no longer works
- `plot(list_of_indicatrices)` has no dispatch

**Plan:**
- Add `"tissot_tbl"` subclass to the tibble returned by `tissot()`, with
  `proj.in` and `proj.out` stored as attributes (not repeated per row)
- Make `indicatrix()` accept either:
  - A `tissot_tbl` object (detect via class; extract proj from attributes)
  - Raw lon/lat + `proj.out` (current signature)
- Move `scale` and `n` to the **plot** method, not `indicatrix()` — these are
  display parameters, not geometric properties
- Return a list with class `"indicatrix_list"` containing:
  - The individual `"indicatrix"` objects (as now)
  - `proj.out` attribute for `tissot_map()` integration

### 4. Restore rich plotting

The old package had detailed multi-element plots:
- Ellipse outline (red fill)
- Lambda direction line (red) — parallel scale
- Phi direction line (blue) — meridional scale
- Major/minor semi-axes
- Reference unit circle (white)

**Current state:** `plot.indicatrix()` only draws a filled polygon.

**Plan:**
- `plot.indicatrix()` — single indicatrix with all axes/reference circle
- `plot.indicatrix_list()` — iterate and plot all, with `add = TRUE` default
- Keep `col`, `border` args; add `show.axes = TRUE`, `show.circle = TRUE`
- Colour-coded fill by distortion metric: `fill.by = "scale.area"` etc.
  using `hcl.colors()` (no extra deps)

### 5. Remove global projection state

`tissot_map()` stores the last projection in `options("tissot.proj.out")`.
Several README examples depend on `tissot_get_proj()`.

**Plan:**
- `tissot_map()` should accept `proj.out` as an explicit argument
- Return it invisibly so it can be captured
- Deprecate (or just remove) `tissot_get_proj()` / `tissot_set_proj()`
- The `indicatrix_list` object carries `proj.out` as an attribute, so
  `plot.indicatrix_list()` can auto-detect the projection for map context

### 6. Set up tests (`testthat`)

**File:** `tests/testthat/test-tissot.R`

Test cases:
- **Conformal projection** (e.g. Stereographic): `angle_deformation ≈ 0`
- **Equal-area projection** (e.g. LAEA): `scale.area ≈ 1`
- **Identity / Plate Carrée at equator**: `scale.h ≈ 1`, `scale.k ≈ 1`
- **Vectorized consistency**: `tissot(c(0,10), c(0,10), ...)` gives same
  results as two separate single-point calls
- **Matrix input**: `tissot(cbind(lon, lat), ...)` same as `tissot(lon, lat, ...)`
- **Indicatrix from tibble**: `indicatrix(tissot(...))` works
- **Indicatrix from lon/lat**: `indicatrix(lon, lat, proj.out=...)` works
- **Ellipse geometry**: `ti_ellipse()` returns closed polygon (first == last row)
- **Edge cases**: Points at poles, antimeridian

Run: `usethis::use_testthat()` then create test files.

### 7. Update README.Rmd

The current README uses several removed/changed APIs:
- `indicatrix0()` — removed, replace with `indicatrix()` single-element access
- `tissot:::.prj()` — removed internal function, replace with
  `gdalraster::transform_xy()`
- `earthcircle::earthcircle()` — not a dependency, remove or replace
- `maps::map()` — not a dependency, use `tissot_map()` or bundled `world` data
- `plot(indicatrix(r, scale=...))` — update to new API
- Several code chunks reference old column names (`axes_*`, `lambda_d*`, `phi_d*`)
  that no longer exist in the new tibble output

**Plan:** Rewrite README to showcase:
1. Basic usage: `tissot()` → tibble of distortion properties
2. Quick plot: `plot(indicatrix(...))` with coastline overlay
3. Interpreting the output: what the columns mean
4. Projection comparison example (2-3 projections side by side)
5. Advanced: accessing individual indicatrix objects, custom plotting

### 8. Clean up man pages

- Remove stale `.Rd` files for `indicatrix0`, `indicatrixes` if they exist
  (check completed — they were already cleaned)
- Expand examples in `tissot.Rd`, `indicatrix.Rd` to show full workflows
- Add `@seealso` cross-references between functions
- Ensure all `@param` docs match current signatures

### 9. Future opportunities (not blocking this refactor)

- `tissot_grid(proj.out, extent, spacing)` — convenience for regular grids
- `summary.tissot_tbl()` — print header with projection info + distortion stats
- `print.tissot_tbl()` — show `"Tissot indicatrix: N points, Robinson"`
- `tissot_graticule()` — draw projected graticule lines
- Natural Earth coastline to replace the `{maps}` derived `world` data
- `geographiclib` (via `geoarea`) in Suggests for geodesic-exact validation
- Colour-coded distortion surface via `image()` / contour
- Projection comparison: `tissot_compare(proj1, proj2, ...)`
- Vignette: `vignette("tissot")` with mathematical background

---

## File inventory (current state on branch)

```
R/tissot.R          — NEW: batched engine, tissot(), indicatrix(), ti_ellipse(), plot.indicatrix()
R/tissot_map.R      — UPDATED: gdalraster::transform_xy, otherwise unchanged
R/tissot-package.R  — UPDATED: import gdalraster
DESCRIPTION         — UPDATED: gdalraster replaces reproj
NAMESPACE           — REGENERATED
man/                — REGENERATED (5 pages: indicatrix, ti_ellipse, tissot, tissot-package, tissot_map)
README.Rmd          — STALE: references old API, needs full rewrite
data/world.rda      — UNCHANGED: bundled coastline matrix
data-raw/DATASET.R  — UNCHANGED
tests/              — DOES NOT EXIST: needs creation
```
