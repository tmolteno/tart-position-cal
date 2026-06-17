"""Core antenna position calibration logic.

Given survey measurements (distances between antennas and to a reference point)
and an initial guess at antenna positions, this module optimizes the positions
to be consistent with the measured distances.
"""

import json
import numpy as np
from scipy.optimize import least_squares, minimize


def geo_angle(x, y):
    """Geographic angle (degrees east of north) from coordinates."""
    return 90.0 - np.degrees(np.arctan2(y, x))


def dist(a, b):
    """Euclidean distance between two points."""
    return np.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2)


def _i_x(i):
    """Index of x-coordinate for antenna i in flat parameter vector."""
    return 2 * i


def _i_y(i):
    """Index of y-coordinate for antenna i in flat parameter vector."""
    return 2 * i + 1


def _p(x, i):
    """Extract (x, y) for antenna i from flat parameter vector."""
    return [x[_i_x(i)], x[_i_y(i)]]


def load_measurements(path):
    """Load antenna distance measurements from an ODS spreadsheet.

    Returns ``(radius, m_ij, n_ant)`` where:
    - *radius*: measured distance from each antenna to the reference point
    - *m_ij*: NxN matrix of inter-antenna distances (NaN where unmeasured)
    - *n_ant*: number of antennas
    """
    import pandas as pd

    cols = [f"A {i}" for i in range(24)]
    data = pd.read_excel(path, "Sheet1", usecols=cols)
    radius = data.loc[0].to_numpy(dtype=float)
    n_ant = len(radius)

    m_ij = np.zeros((n_ant, n_ant))
    for i in range(n_ant):
        m_ij[i, :] = data.loc[i + 1]

    for i in range(n_ant):
        for j in range(n_ant):
            if not np.isnan(m_ij[i, j]):
                m_ij[j, i] = m_ij[i, j]

    return radius, m_ij, n_ant


def load_initial_positions(path):
    """Load initial antenna position guess from a JSON file.

    The JSON must have an ``"antenna_positions"`` key containing a list of
    ``[x, y, z]`` entries in **metres**.  Returns an ``(n_ant, 2)`` array
    in **millimetres** (z is dropped).
    """
    with open(path) as f:
        data = json.load(f)
    positions = np.array(data["antenna_positions"])
    return positions[:, :2] * 1000.0


def fetch_initial_positions(api_url):
    """Fetch initial antenna positions from a TART telescope API.

    Returns an ``(n_ant, 2)`` array in millimetres.
    """
    import requests

    r = requests.get(f"{api_url}/api/v1/imaging/antenna_positions")
    r.raise_for_status()
    positions = np.array(r.json())
    return positions[:, :2] * 1000.0


def calibrate(
    radius,
    m_ij,
    initial_guess,
    *,
    center=(0.0, 0.0),
    rot_index=None,
    rot_degrees=None,
    max_position_error=4200.0,
    radius_weight=10.0,
    rotation_weight=100.0,
    maxiter=500,
    tol=None,
):
    """Run the position calibration optimisation.

    Parameters
    ----------
    radius : (N,) array
        Measured distances from each antenna to the reference point (mm).
    m_ij : (N, N) array
        Matrix of inter-antenna distances (mm); NaN where unmeasured.
    initial_guess : (N, 2) array
        Initial antenna positions in mm.
    center : (float, float)
        Reference point coordinates (mm).
    rot_index : int or None
        Antenna index used to constrain global rotation.  None disables.
    rot_degrees : float or None
        Target geographic angle (degrees) for *rot_index*.
    max_position_error : float
        Half-width of the search bounds around each initial coordinate (mm).
    radius_weight : float
        Weight applied to the radius residual term.
    rotation_weight : float
        Weight applied to the rotation constraint term.
    maxiter : int
        Maximum iterations for the optimiser.
    tol : float or None
        Tolerance for termination.

    Returns
    -------
    result : OptimizeResult
        The scipy optimisation result.  ``result.x`` is the flat parameter
        vector (x0, y0, x1, y1, …) in mm.
    """
    n_ant = len(radius)

    # --- Build the initial flat parameter vector ---
    x0 = np.zeros(2 * n_ant)
    for i in range(n_ant):
        x0[_i_x(i)] = initial_guess[i, 0]
        x0[_i_y(i)] = initial_guess[i, 1]

    # --- Build bounds ---
    bounds = []
    for i in range(n_ant):
        cx, cy = initial_guess[i, 0], initial_guess[i, 1]
        bounds.append((cx - max_position_error, cx + max_position_error))
        bounds.append((cy - max_position_error, cy + max_position_error))

    # --- Pre-compute non-NaN inter-antenna pairs ---
    non_nan_pairs = []
    non_nan_values = []
    for i in range(n_ant):
        for j in range(n_ant):
            if not np.isnan(m_ij[i, j]):
                non_nan_pairs.append((i, j))
                non_nan_values.append(m_ij[i, j])
    non_nan_values = np.array(non_nan_values)

    # --- Residual functions ---
    def radius_residual(x):
        pred = np.array([dist(center, _p(x, i)) for i in range(n_ant)])
        return pred - radius

    def m_ij_residual(x):
        pred = np.array(
            [dist(_p(x, i), _p(x, j)) for (i, j) in non_nan_pairs]
        )
        return pred - non_nan_values

    def rot_residual(x):
        if rot_index is None or rot_degrees is None:
            return 0.0
        xi, yi = _p(x, rot_index)
        return geo_angle(xi, yi) - rot_degrees

    def objective(x):
        val = np.sum(radius_residual(x) ** 2) * radius_weight
        val += np.sum(m_ij_residual(x) ** 2)
        val += (rot_residual(x) ** 2) * rotation_weight
        return val

    options = {"maxiter": maxiter}
    if tol is not None:
        options["xatol"] = tol
        options["fatol"] = tol

    return minimize(objective, x0, bounds=bounds, options=options)


def _tukey_weights(residuals, scale=4.685):
    """Tukey's bisquare weights for IRLS.

    *scale* is the tuning constant (default 4.685 for 95% efficiency
    under normality).  Residuals are Studentised internally via MAD.
    """
    mad = np.median(np.abs(residuals - np.median(residuals)))
    if mad < 1e-12:
        return np.ones_like(residuals)
    u = residuals / (scale * mad * 1.4826)
    w = np.where(np.abs(u) < 1.0, (1.0 - u**2) ** 2, 0.0)
    return np.clip(w, 0.0, 1.0)


def _huber_weights(residuals, scale=1.345):
    """Huber weights for IRLS.

    *scale* is the tuning constant (default 1.345 for 95% efficiency
    under normality).
    """
    mad = np.median(np.abs(residuals - np.median(residuals)))
    if mad < 1e-12:
        return np.ones_like(residuals)
    u = np.abs(residuals) / (scale * mad * 1.4826)
    w = np.where(u <= 1.0, 1.0, 1.0 / u)
    return np.clip(w, 0.0, 1.0)


def calibrate_irls(
    radius,
    m_ij,
    initial_guess,
    *,
    center=(0.0, 0.0),
    rot_index=None,
    rot_degrees=None,
    max_position_error=4200.0,
    radius_weight=10.0,
    rotation_weight=100.0,
    maxiter=500,
    irls_max_iter=10,
    irls_tol=1e-4,
    weight_function="tukey",
):
    """Run calibration using Iteratively Reweighted Least Squares (IRLS).

    IRLS iteratively solves a weighted least-squares problem, recomputing
    per-measurement weights at each iteration to down-weight outliers.
    This is more robust than the standard ``calibrate`` when the survey
    data contains erroneous measurements.

    Parameters are the same as :func:`calibrate` with these additions:

    irls_max_iter : int
        Maximum number of IRLS outer iterations (default 10).
    irls_tol : float
        Convergence tolerance on the parameter vector between IRLS
        iterations.
    weight_function : str
        ``"tukey"`` (Tukey's bisquare) or ``"huber"``.
    """
    if weight_function == "tukey":
        weight_fn = _tukey_weights
    elif weight_function == "huber":
        weight_fn = _huber_weights
    else:
        raise ValueError(f"Unknown weight_function: {weight_function}")

    n_ant = len(radius)

    # --- Build initial flat parameter vector ---
    x0 = np.zeros(2 * n_ant)
    for i in range(n_ant):
        x0[_i_x(i)] = initial_guess[i, 0]
        x0[_i_y(i)] = initial_guess[i, 1]

    # --- Build bounds ---
    bounds_lower = []
    bounds_upper = []
    for i in range(n_ant):
        cx, cy = initial_guess[i, 0], initial_guess[i, 1]
        bounds_lower.append(cx - max_position_error)
        bounds_lower.append(cy - max_position_error)
        bounds_upper.append(cx + max_position_error)
        bounds_upper.append(cy + max_position_error)
    bounds = (bounds_lower, bounds_upper)

    # --- Pre-compute non-NaN inter-antenna pairs ---
    non_nan_pairs = []
    non_nan_values = []
    for i in range(n_ant):
        for j in range(n_ant):
            if not np.isnan(m_ij[i, j]):
                non_nan_pairs.append((i, j))
                non_nan_values.append(m_ij[i, j])
    non_nan_values = np.array(non_nan_values)

    n_radius = n_ant
    n_ij = len(non_nan_values)
    n_rot = 1 if (rot_index is not None and rot_degrees is not None) else 0
    n_total = n_radius + n_ij + n_rot

    sqrt_radius_weight = np.sqrt(radius_weight)
    sqrt_rotation_weight = np.sqrt(rotation_weight)

    def combined_residuals(x):
        """Flat residual vector (unweighted)."""
        res = np.zeros(n_total)
        # Radius residuals
        for i in range(n_ant):
            res[i] = sqrt_radius_weight * (dist(center, _p(x, i)) - radius[i])
        # Inter-antenna residuals
        for k, (i, j) in enumerate(non_nan_pairs):
            res[n_radius + k] = dist(_p(x, i), _p(x, j)) - non_nan_values[k]
        # Rotation residual
        if n_rot:
            xi, yi = _p(x, rot_index)
            res[-1] = sqrt_rotation_weight * (geo_angle(xi, yi) - rot_degrees)
        return res

    # --- IRLS outer loop ---
    weights = np.ones(n_total)
    prev_x = x0.copy()

    for irls_iter in range(irls_max_iter):
        def weighted_residuals(x):
            return np.sqrt(weights) * combined_residuals(x)

        res_ls = least_squares(
            weighted_residuals,
            prev_x,
            bounds=bounds,
            method="trf",
            max_nfev=maxiter,
            xtol=1e-12,
            ftol=1e-12,
        )
        new_x = res_ls.x

        # Update weights from the *unweighted* residuals
        raw_res = combined_residuals(new_x)
        weights = weight_fn(raw_res)

        shift = np.max(np.abs(new_x - prev_x))
        prev_x = new_x

        if shift < irls_tol:
            break

    # Build a MinimizeResult-compatible result
    res_ls.x = prev_x
    res_ls.fun = np.sum(combined_residuals(prev_x) ** 2)
    res_ls.nfev = getattr(res_ls, "nfev", 0)
    res_ls.njev = getattr(res_ls, "njev", 0)
    res_ls.nit = getattr(res_ls, "nit", 0)
    return res_ls


def result_to_positions(result, n_ant):
    """Extract (n_ant, 2) positions array from an OptimizeResult."""
    return result.x.reshape((n_ant, 2))


def result_to_json(result, n_ant):
    """Convert an OptimizeResult to the API-compatible JSON dict.

    Positions are rounded to 3 decimal places (metres) with z=0.
    """
    pos_mm = result_to_positions(result, n_ant)
    out = np.zeros((n_ant, 3))
    out[:, :2] = np.round(pos_mm / 1000.0, 3)
    return {"antenna_positions": out.tolist()}


def compute_residuals(result, radius, m_ij, n_ant):
    """Return (radius_residuals, ij_residuals, ij_pairs) for diagnostics."""
    x = result.x

    r_res = np.array(
        [dist((0, 0), _p(x, i)) - radius[i] for i in range(n_ant)]
    )

    pairs = []
    ij_res = []
    for i in range(n_ant):
        for j in range(n_ant):
            if not np.isnan(m_ij[i, j]):
                pairs.append((i, j))
                ij_res.append(dist(_p(x, i), _p(x, j)) - m_ij[i, j])

    return r_res, np.array(ij_res), np.array(pairs)


def plot_positions(initial, final, outpath):
    """Save a side-by-side plot of initial vs final positions."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))

    ax[0].scatter(initial[:, 0], initial[:, 1], label="initial", color="blue")
    for i in range(initial.shape[0]):
        ax[0].text(initial[i, 0], initial[i, 1], str(i), color="orange")
    ax[0].set_title("Initial Guess")
    ax[0].grid(True)
    ax[0].set_xlabel("x (mm)")
    ax[0].set_ylabel("y (mm)")

    ax[1].scatter(final[:, 0], final[:, 1], label="final", color="orange")
    for i in range(final.shape[0]):
        ax[1].text(final[i, 0], final[i, 1], str(i), color="blue")
    ax[1].set_title("Final Solution")
    ax[1].grid(True)
    ax[1].set_xlabel("x (mm)")
    ax[1].set_ylabel("y (mm)")

    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)


def plot_differences(initial, final, outpath):
    """Save a plot of position differences (final - initial)."""
    import matplotlib.pyplot as plt

    diff = final - initial[:, :2]

    fig, ax = plt.subplots(figsize=(12, 12))
    ax.scatter(diff[:, 0], diff[:, 1], color="red")
    for i in range(diff.shape[0]):
        ax.text(diff[i, 0], diff[i, 1], str(i))
    ax.grid(True)
    ax.set_title("Differences from initial position")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    fig.savefig(outpath)
    plt.close(fig)


def plot_residual_histogram(ij_residuals, outpath):
    """Save a histogram of inter-antenna distance residuals."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    ax.hist(ij_residuals, bins="fd")
    ax.set_title("Histogram of residuals")
    ax.set_xlabel("Residual (mm)")
    ax.grid(True)
    fig.savefig(outpath)
    plt.close(fig)


def print_residual_report(result, radius, m_ij, n_ant):
    """Print a residual summary to stdout."""
    r_res, ij_res, pairs = compute_residuals(result, radius, m_ij, n_ant)

    print("Radius residuals (mm):")
    for i in range(n_ant):
        print(f"  Ant {i:2d}: {r_res[i]:+7.2f}")

    pct = 90
    p95 = np.percentile(np.abs(ij_res), pct)
    print(f"\n{pct}th percentile of |ij residuals|: {p95:.2f} mm")

    big = []
    for r, (i, j) in zip(ij_res, pairs):
        if np.abs(r) > p95 and i > j:
            big.append((r, i, j))

    if big:
        print(f"\nLargest residuals (>p{pct}):")
        for r, i, j in sorted(big):
            print(f"  res[{i},{j}] = {r:+.1f}")
