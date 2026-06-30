# Changelog

## Unreleased

## [0.2.0] - 2026-06-30

### Added

- `compute_bearing(lat1, lon1, lat2, lon2)` function: solves the WGS84
  inverse geodesic problem via PROJ (pyproj) and returns the geographic
  bearing (degrees east of north) from the phase-centre to a distant
  landmark.  Exported from `tart_position_cal`.
- `--chirality-index` CLI argument: breaks the reflection (chirality)
  degeneracy by enforcing the sign of the cross product from the initial
  position estimates.  Documented in `doc/SYMMETRY.md`.
- `--irls-weight-fn`, `--irls-max-iter`, `--plot-dir`, `--title`
  documented in paper and README.
- `doc/SYMMETRY.md`: documents the three symmetries of the distance-data
  problem (translation, rotation, reflection) and how each is resolved.
- `doc/tart_position_cal.tex` §2.4: symmetries and global frame
  determination section in the journal article.
- Input validation in `calibrate` / `calibrate_irls` (via the shared
  `_build_problem`): raises `ValueError` with a clear message when
  `initial_guess` or `m_ij` dimensions disagree with `len(radius)`, or when
  `rot_index` / `chirality_index` are outside `[0, n_ant)`. Previously these
  surfaced as opaque `IndexError`s deep in the optimiser.
- `fetch_initial_positions(api_url, *, timeout=30.0)` now passes a network
  timeout to `requests.get`, so an unresponsive API fails fast instead of
  hanging indefinitely.  Exposed on the CLI as `--api-timeout` (default 30 s).
- `load_measurements` / `load_initial_positions` now raise clear, actionable
  errors for bad inputs: a missing file gives a `FileNotFoundError` naming
  the path, malformed JSON and missing `antenna_positions` keys give
  `ValueError`s, and a measurements spreadsheet without a `Sheet1` sheet
  (or too few rows for the antenna count) reports what was expected instead
  of a raw pandas/`KeyError` traceback.
- CLI: moved the previously buried `import os` (formerly inline inside
  `main()`, after its first use) to the module imports at the top of
  `tart_position_cal.py`.
- `tests/test_calibrate.py`: 31 tests covering `compute_bearing`,
  `geo_angle`, `calibrate`, `calibrate_irls`, chirality constraint,
  bounds, collinearity guard, input validation, the API timeout, and loader
  error handling.
- Tool configuration in `pyproject.toml`: `[tool.ruff]` pins
  `line-length = 100` and `target-version = "py311"`, with `[tool.ruff.lint]`
  `extend-select = ["E501"]` enforcing the line length on top of ruff's
  defaults; `[tool.pytest.ini_options]` sets `testpaths = ["tests"]`.

### Changed

- `doc/tart_position_cal.tex` synced to the 0.2.0 algorithm: the paper now
  describes the trust-region reflective solver (not L-BFGS-B), and the
  objective (Eq.~2) shows the radius weight $w_r$, the unique-pair sum
  ($\sum_{i<j}$, 276 for $N=24$), and the chirality term $L_\chi$.  All
  UNAM results (\S4) were regenerated from the current code: final objective
  $4.35\times10^2$ (20 iterations), residual MAD/std/p50/p90 of
  0.63/1.11/0.59/1.65~mm, the IRLS comparison (0.49/1.11/0.49/1.44~mm in 4
  iterations, with the unweighted objective rising to $5.40\times10^2$), and
  the largest-residuals table.
- Lint cleanup: `ruff check .` now passes cleanly. Removed an unused
  `import math`, consolidated three scattered `from tart_position_cal...`
  imports in the test module into one top-level import (fixing E402), and
  dropped stray `f` prefixes from two format strings without placeholders
  (fixing F541).
- `calibrate` and `calibrate_irls` now share a single `_Problem` setup
  (initial vector, bounds, pair enumeration, chirality bookkeeping) and a
  common residual-vector evaluator, removing the duplication between them.
- The standard `calibrate` solver switched from `scipy.optimize.minimize`
  (L-BFGS-B, scalar objective) to `scipy.optimize.least_squares`
  (method='trf'), matching `calibrate_irls`.  This exploits the residual
  Jacobian and converges in far fewer evaluations (~23 vs ~3920 fevals on
  the UNAM example).  `result.fun` is still the scalar sum-of-squares
  objective; `result.nit` now holds the Jacobian-evaluation count.
- Inter-antenna distance pairs are now counted once each (upper triangle)
  instead of twice, in the objective and in the diagnostic report.  This
  removes an unintended 2× weighting of the inter-antenna term relative to
  the radius term.  Calibrated positions shift by ≤ 1 mm on the UNAM data.
- IRLS mode: rotation and chirality constraints are now excluded from
  Tukey/Huber reweighting (they are priors, not noisy measurements).
- Equation (3) in paper: `arctan(x/y)` replaced with `arctan2(x, y)`.
- Flow-chart `overlay` option removed to prevent clipping artifacts.

### Fixed

- `calibrate(tol=...)` was silently ignored: it passed `xatol`/`fatol`, which
  L-BFGS-B does not recognise (`OptimizeWarning: Unknown solver options`).
  `tol` now sets `ftol`/`xtol`/`gtol` of the trust-region least-squares
  solver and actually controls convergence.
- Garbled equations (`\delt^sa_`, `L_\theta(...)\theta_k`) from Kirsten's
  grammar restructure.
- Typo: "atendees" → "attendees", "fewer measurement" → "fewer measurements".
- Orphaned "Where" and disconnected sentence about θ_k in §2.1.
- Removed "TKTK" placeholder comment.
- Added `\label{fig:spreadsheet}` and `\label{sec:irls}` for cross-references.
