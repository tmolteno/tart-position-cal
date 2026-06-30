# Changelog

## Unreleased

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
- `tests/test_calibrate.py`: 16 tests covering `compute_bearing`,
  `geo_angle`, `calibrate`, `calibrate_irls`, chirality constraint,
  bounds, and collinearity guard.

### Changed

- IRLS mode: rotation and chirality constraints are now excluded from
  Tukey/Huber reweighting (they are priors, not noisy measurements).
- `scipy.optimize.minimize` replaced with `scipy.optimize.least_squares`
  (method='trf') for the IRLS inner solver.
- Equation (3) in paper: `arctan(x/y)` replaced with `arctan2(x, y)`.
- Flow-chart `overlay` option removed to prevent clipping artifacts.

### Fixed

- Garbled equations (`\delt^sa_`, `L_\theta(...)\theta_k`) from Kirsten's
  grammar restructure.
- Typo: "atendees" → "attendees", "fewer measurement" → "fewer measurements".
- Orphaned "Where" and disconnected sentence about θ_k in §2.1.
- Removed "TKTK" placeholder comment.
- Added `\label{fig:spreadsheet}` and `\label{sec:irls}` for cross-references.
