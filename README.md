# tart-position-cal

Calibrate TART radio telescope antenna positions from survey measurements.

## Installation

```bash
uv tool install tart-position-cal
```

Or from source:

```bash
cd tart-position-cal
uv sync
```

## Usage

### Obtaining the reference bearing

The geographic bearing of the reference antenna is computed from the
phase-centre and landmark coordinates using `compute_bearing()`:

```python
from tart_position_cal import compute_bearing

# UNAM: phase-centre to distant hill
bearing = compute_bearing(-22.612053, 17.056784,
                          -22.553822, 17.077290)
# bearing = 18.112521...
```

### Running the calibration

```bash
tart-position-cal \
    --measurements antenna_measurements.ods \
    --positions antenna_positions.json \
    --output calibrated_antenna_positions.json \
    --rot-index 3 \
    --rot-degrees 18.1125 \
    --chirality-index 10
```

### Options

| Option | Description |
|--------|-------------|
| `--measurements` | Path to .ods file with antenna distance measurements (required) |
| `--positions` | Path to JSON file with initial antenna positions in metres |
| `--api` | TART API base URL to fetch initial positions from |
| `--api-timeout` | Network timeout in seconds for `--api` requests (default: 30) |
| `--output` | Path for output calibrated positions JSON (required) |
| `--rot-index` | Antenna index for global rotation constraint (default: 0) |
| `--rot-degrees` | Target geographic angle in degrees for constrained antenna (default: 0) |
| `--chirality-index` | Antenna index to break reflection (chirality) degeneracy. Requires `--rot-index`. |
| `--max-iter` | Maximum optimizer iterations (default: 500) |
| `--max-position-error` | Search bound half-width in mm (default: 4200) |
| `--no-plots` | Skip diagnostic plot generation |
| `--plot-dir` | Directory for diagnostic plots (default: current dir) |
| `--irls` | Use Iteratively Reweighted Least Squares for robust calibration |
| `--irls-weight-fn` | Weight function for IRLS: `tukey` (bisquare) or `huber` (default: `tukey`) |
| `--irls-max-iter` | Maximum IRLS outer iterations (default: 10) |
| `--title` | Optional prefix prepended to all output file names (e.g. `na-unam`) |

### Input format

**Measurements** (.ods spreadsheet):
- Row 0: distance from each antenna to the reference point
- Rows 1..N: inter-antenna distance matrix

**Initial positions** (JSON):
```json
{
    "antenna_positions": [
        [x0, y0, z0],
        [x1, y1, z1],
        ...
    ]
}
```

Positions are in metres; the calibration works in millimetres internally.

### Output

A JSON file with calibrated 3D positions (z=0) in metres, suitable for upload to the TART telescope API:

```bash
tart_upload_antenna_positions --api https://api.elec.ac.nz/tart/site-name --file calibrated_antenna_positions.json --pw XXXX
```

### Robust calibration with IRLS

When survey measurements contain outliers (erroneous distance readings), standard
least-squares optimisation can be skewed by those bad measurements.  Pass `--irls`
to use Iteratively Reweighted Least Squares:

```bash
tart-position-cal \
    --measurements antenna_measurements.ods \
    --positions antenna_positions.json \
    --output calibrated_antenna_positions.json \
    --rot-index 3 \
    --rot-degrees 18.1125 \
    --chirality-index 10 \
    --irls
```

IRLS solves a series of weighted least-squares problems, recomputing per-measurement
weights after each solve to down-weight outliers.  Two weight functions are available:

| `--irls-weight-fn` | Behaviour |
|--------------------|-----------|
| `tukey` (default)  | Tukey's bisquare — zeros out clear outliers entirely |
| `huber`            | Huber — caps the influence of outliers with a saturating penalty |

The outer loop stops when the parameter vector stops moving (`< 1e-4` mm) or after
`--irls-max-iter` iterations (default 10).

## Example: UNAM site (Namibia)

The TART station at UNAM, Namibia (24 antennas on 5 spiral arms).

**Bearing:** antenna 3 aligns with a distant hill.  Compute the bearing:

```python
from tart_position_cal import compute_bearing
bearing = compute_bearing(-22.612053, 17.056784,
                          -22.553822, 17.077290)
# bearing = 18.112521...
```

**Calibration:**

```bash
tart-position-cal \
    --measurements example/antenna_measurements.ods \
    --positions example/antenna_positions.json \
    --output na-unam_processed_antenna_positions.json \
    --rot-index 3 \
    --rot-degrees 18.1125 \
    --chirality-index 10
```

```
Loading measurements from example/antenna_measurements.ods
  Found 24 antennas
Loading initial positions from example/antenna_positions.json
Running calibration (rot_index=3, rot_degrees=18.1125)
  Optimisation finished: `ftol` termination condition is satisfied.
  Final objective: 434.6330
  Iterations: 20
  Function evaluations: 23
Calibrated positions written to na-unam_processed_antenna_positions.json
Radius residuals (mm):
  Ant  0:   +0.97
  Ant  1:   +0.63
  Ant  2:   +1.16
  Ant  3:   +0.31
  Ant  4:   +0.86
  Ant  5:   -0.07
  Ant  6:   +0.22
  Ant  7:   -0.06
  Ant  8:   +0.22
  Ant  9:   -0.18
  Ant 10:   +0.48
  Ant 11:   -0.21
  Ant 12:   +0.90
  Ant 13:   +1.40
  Ant 14:   +0.65
  Ant 15:   -0.12
  Ant 16:   +0.39
  Ant 17:   +0.23
  Ant 18:   +0.67
  Ant 19:   +0.75
  Ant 20:   +0.08
  Ant 21:   -0.10
  Ant 22:   +0.33
  Ant 23:   -0.02

Inter-antenna residual statistics:
  Median Absolute Deviation: 0.63 mm
  Standard Deviation:       1.11 mm
  50th percentile |res|:    0.59 mm
  90th percentile |res|:    1.65 mm

Largest residuals (>90th pct):
  res[2,13] = -4.0
  res[3,7] = -3.5
  res[1,14] = -3.5
  res[4,19] = -3.1
  res[8,17] = -2.8
  res[0,12] = -2.5
  res[2,5] = -2.5
  res[4,14] = -2.5
  res[0,14] = -2.4
  res[3,12] = -2.3
  res[9,22] = -2.3
  res[13,22] = -2.2
  res[0,13] = -2.1
  res[9,17] = -1.9
  res[2,12] = -1.9
  res[7,13] = -1.8
  res[1,13] = -1.8
  res[9,15] = -1.8
  res[1,12] = -1.7
  res[16,20] = -1.7
  res[3,15] = +1.8
  res[15,23] = +1.9
  res[15,22] = +2.2
  res[12,22] = +3.1
  res[4,17] = +3.4
  res[2,9] = +4.1
  res[16,17] = +4.5
  res[9,19] = +5.0
```

The 90th-percentile inter-antenna residual of 1.65 mm shows the survey
measurements are consistent to about ±1.7 mm.  The largest residuals (up to
~5 mm) flag a few measurements that may warrant re-surveying.

## Development

```bash
uv sync --group dev
uv run pytest
```
