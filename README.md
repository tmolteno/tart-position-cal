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

```bash
tart-position-cal \
    --measurements antenna_measurements.ods \
    --positions antenna_positions.json \
    --output calibrated_antenna_positions.json \
    --rot-index 3 \
    --rot-degrees 18.1125
```

### Options

| Option | Description |
|--------|-------------|
| `--measurements` | Path to .ods file with antenna distance measurements (required) |
| `--positions` | Path to JSON file with initial antenna positions in metres |
| `--api` | TART API base URL to fetch initial positions from |
| `--output` | Path for output calibrated positions JSON (required) |
| `--rot-index` | Antenna index for global rotation constraint (default: 0) |
| `--rot-degrees` | Target geographic angle in degrees for constrained antenna (default: 0) |
| `--max-iter` | Maximum optimizer iterations (default: 500) |
| `--max-position-error` | Search bound half-width in mm (default: 4200) |
| `--no-plots` | Skip diagnostic plot generation |
| `--plot-dir` | Directory for diagnostic plots (default: current dir) |
| `--irls` | Use Iteratively Reweighted Least Squares for robust calibration |
| `--irls-weight-fn` | Weight function for IRLS: `tukey` (bisquare) or `huber` (default: `tukey`) |
| `--irls-max-iter` | Maximum IRLS outer iterations (default: 10) |

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

## Example: na-unam site

The TART station at UNAM, Namibia (24 antennas on 5 spiral arms).  A rotation
constraint fixes antenna 3 at 18.1125° geographic to resolve the global orientation
ambiguity.

```bash
tart-position-cal \
    --measurements example/antenna_measurements.ods \
    --positions example/antenna_positions.json \
    --output na-unam_processed_antenna_positions.json \
    --rot-index 3 \
    --rot-degrees 18.1125
```

```
Loading measurements from example/antenna_measurements.ods
  Found 24 antennas
Loading initial positions from example/antenna_positions.json
Running calibration (rot_index=3, rot_degrees=18.1125)
  Optimisation finished: CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH
  Final objective: 765.8097
  Iterations: 71
  Function evaluations: 3920
Calibrated positions written to na-unam_processed_antenna_positions.json
Radius residuals (mm):
  Ant  0:   +1.24
  Ant  1:   +0.79
  Ant  2:   +1.43
  Ant  3:   +0.36
  Ant  4:   +1.03
  Ant  5:   -0.12
  Ant  6:   +0.24
  Ant  7:   -0.12
  Ant  8:   +0.24
  Ant  9:   -0.25
  Ant 10:   +0.61
  Ant 11:   -0.28
  Ant 12:   +1.10
  Ant 13:   +1.70
  Ant 14:   +0.78
  Ant 15:   -0.17
  Ant 16:   +0.47
  Ant 17:   +0.26
  Ant 18:   +0.80
  Ant 19:   +0.89
  Ant 20:   +0.06
  Ant 21:   -0.17
  Ant 22:   +0.36
  Ant 23:   -0.06

90th percentile of |ij residuals|: 1.53 mm

Largest residuals (>p90):
  res[7,3] = -3.5
  res[13,2] = -3.4
  res[14,1] = -3.2
  res[17,8] = -2.8
  res[19,4] = -2.8
  res[22,9] = -2.3
  res[5,2] = -2.2
  res[14,4] = -2.2
  res[12,3] = -2.1
  res[12,0] = -2.0
  res[14,0] = -2.0
  res[17,9] = -2.0
  res[15,9] = -1.9
  res[22,13] = -1.9
  res[13,7] = -1.7
  res[20,16] = -1.6
  res[13,0] = -1.5
  res[22,21] = -1.5
  res[11,7] = +1.6
  res[22,16] = +1.6
  res[23,15] = +1.8
  res[15,3] = +1.8
  res[22,15] = +2.2
  res[22,12] = +3.3
  res[17,4] = +3.6
  res[9,2] = +4.2
  res[17,16] = +4.4
  res[19,9] = +5.0
```

The 90th-percentile inter-antenna residual of 1.53 mm shows the survey
measurements are consistent to about ±1.5 mm.  The largest residuals (up to
~5 mm) flag a few measurements that may warrant re-surveying.

## Development

```bash
uv sync --group dev
uv run pytest
```
