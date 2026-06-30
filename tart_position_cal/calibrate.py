"""Core antenna position calibration logic.

Given survey measurements (distances between antennas and to a reference point)
and an initial guess at antenna positions, this module optimizes the positions
to be consistent with the measured distances.
"""

import json
from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares


def geo_angle(x, y):
    """Geographic angle (degrees east of north) from coordinates."""
    return 90.0 - np.degrees(np.arctan2(y, x))


def dist(a, b):
    """Euclidean distance between two points."""
    return np.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2)


def compute_bearing(lat1, lon1, lat2, lon2):
    """Compute the geographic bearing from one point to another.

    Solves the inverse geodesic problem on the WGS84 ellipsoid using the
    PROJ library.  Returns the forward azimuth (degrees east of north)
    from ``(lat1, lon1)`` to ``(lat2, lon2)``.

    This is the value to pass to ``--rot-degrees``.

    Parameters
    ----------
    lat1, lon1 : float
        Latitude and longitude of the origin (the TART phase-centre),
        in decimal degrees.
    lat2, lon2 : float
        Latitude and longitude of the destination (the distant landmark),
        in decimal degrees.

    Returns
    -------
    bearing : float
        Geographic bearing in degrees east of north.

    Examples
    --------
    >>> # UNAM telescope to distant hill
    >>> compute_bearing(-22.612053, 17.056784, -22.553822, 17.077290)
    18.1125...
    """
    try:
        import pyproj
    except ImportError:
        raise ImportError(
            "compute_bearing requires pyproj.  Install it with: pip install pyproj"
        )
    geodesic = pyproj.Geod(ellps="WGS84")
    fwd_azimuth, _back_azimuth, _distance = geodesic.inv(lon1, lat1, lon2, lat2)
    return fwd_azimuth


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
    try:
        data = pd.read_excel(path, "Sheet1", usecols=cols)
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            f"Cannot read measurements from {path!r}: file not found"
        ) from exc
    except ValueError as exc:
        # pandas raises ValueError for a missing sheet or an unreadable format.
        raise ValueError(
            f"Could not read measurements from {path!r}: {exc}. "
            f"Expected an ODS spreadsheet with a sheet named 'Sheet1'."
        ) from exc

    if len(data) < 1:
        raise ValueError(f"Measurements file {path!r} has no data rows")
    radius = data.loc[0].to_numpy(dtype=float)
    n_ant = len(radius)

    if len(data) < n_ant + 1:
        raise ValueError(
            f"Measurements file {path!r} has {len(data)} data row(s) but the "
            f"radius row implies {n_ant} antennas; expected at least "
            f"{n_ant + 1} rows (1 radius + {n_ant} inter-antenna distances)."
        )

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
    try:
        with open(path) as f:
            data = json.load(f)
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            f"Cannot read initial positions from {path!r}: file not found"
        ) from exc
    except json.JSONDecodeError as exc:
        raise ValueError(
            f"Invalid JSON in initial-positions file {path!r}: "
            f"{exc.msg} (line {exc.lineno}, column {exc.colno})"
        ) from exc

    if "antenna_positions" not in data:
        raise ValueError(
            f"Initial-positions file {path!r} has no 'antenna_positions' key"
        )
    positions = np.array(data["antenna_positions"])
    return positions[:, :2] * 1000.0


def fetch_initial_positions(api_url, *, timeout=30.0):
    """Fetch initial antenna positions from a TART telescope API.

    Parameters
    ----------
    api_url : str
        Base URL of the TART API.
    timeout : float or tuple
        Network timeout in seconds passed to :func:`requests.get`.  May be a
        single value or a ``(connect, read)`` tuple.  Defaults to 30 s so an
        unresponsive API fails fast instead of hanging indefinitely.

    Returns an ``(n_ant, 2)`` array in millimetres.
    """
    import requests

    r = requests.get(
        f"{api_url}/api/v1/imaging/antenna_positions",
        timeout=timeout,
    )
    r.raise_for_status()
    positions = np.array(r.json())
    return positions[:, :2] * 1000.0


# ---------------------------------------------------------------------------
# Shared problem setup (used by both calibrate and calibrate_irls)
# ---------------------------------------------------------------------------


@dataclass
class _Problem:
    """Precomputed calibration problem shared by the solvers.

    Bundles the initial parameter vector, bounds, the unique inter-antenna
    distance pairs, and the rotation/chirality constraint bookkeeping, and
    evaluates the flat residual vector consumed by the least-squares solvers.
    """

    n_ant: int
    x0: np.ndarray
    lower: np.ndarray
    upper: np.ndarray
    center: tuple[float, float]
    radius: np.ndarray
    pairs: list[tuple[int, int]]
    pair_values: np.ndarray
    sqrt_radius_weight: float
    sqrt_rotation_weight: float
    sqrt_chirality_weight: float
    rot_index: object = None
    rot_degrees: object = None
    chirality_index: object = None
    expected_chirality_sign: object = None
    chirality_norm: object = None

    @property
    def bounds(self):
        return (self.lower, self.upper)

    @property
    def n_radius(self):
        return self.n_ant

    @property
    def n_ij(self):
        return len(self.pair_values)

    @property
    def n_rot(self):
        return 1 if (self.rot_index is not None and self.rot_degrees is not None) else 0

    @property
    def n_chirality(self):
        return (
            1
            if (self.chirality_index is not None and self.rot_index is not None)
            else 0
        )

    @property
    def n_total(self):
        return self.n_radius + self.n_ij + self.n_rot + self.n_chirality

    def combined_residuals(self, x):
        """Flat residual vector with term weights baked in (unweighted by IRLS)."""
        res = np.zeros(self.n_total)
        # Radius residuals
        for i in range(self.n_ant):
            res[i] = self.sqrt_radius_weight * (
                dist(self.center, _p(x, i)) - self.radius[i]
            )
        # Inter-antenna residuals (one entry per unique pair)
        for k, (i, j) in enumerate(self.pairs):
            res[self.n_radius + k] = dist(_p(x, i), _p(x, j)) - self.pair_values[k]
        # Rotation residual
        if self.n_rot:
            xi, yi = _p(x, self.rot_index)
            res[self.n_radius + self.n_ij] = self.sqrt_rotation_weight * (
                geo_angle(xi, yi) - self.rot_degrees
            )
        # Chirality residual
        if self.n_chirality:
            x_ref, y_ref = _p(x, self.rot_index)
            x_chi, y_chi = _p(x, self.chirality_index)
            cross = x_ref * y_chi - y_ref * x_chi
            normalized = cross / self.chirality_norm
            res[self.n_radius + self.n_ij + self.n_rot] = (
                self.sqrt_chirality_weight
                * min(0.0, normalized * self.expected_chirality_sign)
            )
        return res


def _build_problem(
    radius,
    m_ij,
    initial_guess,
    *,
    center=(0.0, 0.0),
    rot_index=None,
    rot_degrees=None,
    chirality_index=None,
    chirality_sign=None,
    chirality_weight=1e6,
    max_position_error=4200.0,
    radius_weight=10.0,
    rotation_weight=100.0,
):
    """Pre-compute everything needed to evaluate the calibration residuals.

    Returns a :class:`_Problem` carrying the initial parameter vector,
    bounds, the unique inter-antenna pairs, and the rotation/chirality
    constraint bookkeeping.  Both :func:`calibrate` and :func:`calibrate_irls`
    build on top of the same setup.
    """
    n_ant = len(radius)

    # --- Validate inputs (fail early with clear messages) ---
    initial_guess = np.asarray(initial_guess, dtype=float)
    if initial_guess.ndim != 2 or initial_guess.shape[1] != 2:
        raise ValueError(
            f"initial_guess must have shape (n_ant, 2); got {initial_guess.shape}"
        )
    if initial_guess.shape[0] != n_ant:
        raise ValueError(
            f"initial_guess has {initial_guess.shape[0]} rows but len(radius)={n_ant}"
        )
    m_ij = np.asarray(m_ij, dtype=float)
    if m_ij.shape != (n_ant, n_ant):
        raise ValueError(f"m_ij must have shape ({n_ant}, {n_ant}); got {m_ij.shape}")
    for _name, _idx in (("rot_index", rot_index), ("chirality_index", chirality_index)):
        if _idx is not None and not (0 <= _idx < n_ant):
            raise ValueError(
                f"{_name}={_idx} is out of range for {n_ant} antennas "
                f"(must be 0..{n_ant - 1})"
            )

    # --- Initial flat parameter vector: (x0, y0, x1, y1, ...) ---
    x0 = np.empty(2 * n_ant)
    x0[0::2] = initial_guess[:, 0]
    x0[1::2] = initial_guess[:, 1]

    # --- Bounds: box of half-width max_position_error around each antenna ---
    lower = np.empty(2 * n_ant)
    upper = np.empty(2 * n_ant)
    lower[0::2] = initial_guess[:, 0] - max_position_error
    lower[1::2] = initial_guess[:, 1] - max_position_error
    upper[0::2] = initial_guess[:, 0] + max_position_error
    upper[1::2] = initial_guess[:, 1] + max_position_error

    # --- Unique inter-antenna pairs (upper triangle, excluding diagonal) ---
    # load_measurements symmetrises m_ij, so each measured pair appears at
    # both (i, j) and (j, i).  Taking i < j keeps each pair once, avoiding
    # double-counting it in the objective and diagnostics.
    pairs = [
        (i, j)
        for i in range(n_ant)
        for j in range(i + 1, n_ant)
        if not np.isnan(m_ij[i, j])
    ]
    pair_values = np.array([m_ij[i, j] for (i, j) in pairs])

    # --- Chirality sign: declared explicitly by the caller, NOT inferred ---
    # from the initial guess.  The initial guess is only used here to check
    # that (rot_index, chirality_index) is not a degenerate (collinear) pair;
    # it is deliberately *not* used to decide which sign to enforce, because
    # an as-built array can be mirrored relative to its nominal design file
    # (this happened at UNAM) and the optimiser cannot detect that on its
    # own.  The caller must determine the true sign from the physical array
    # (see doc/SYMMETRY.md) and pass it in.
    expected_chirality_sign = None
    chirality_norm = None
    if chirality_index is not None and rot_index is not None:
        x_ref = initial_guess[rot_index, 0]
        y_ref = initial_guess[rot_index, 1]
        x_chi = initial_guess[chirality_index, 0]
        y_chi = initial_guess[chirality_index, 1]
        cross = x_ref * y_chi - y_ref * x_chi
        if abs(cross) < 1e-6:
            raise ValueError(
                f"chirality_index={chirality_index} is collinear with "
                f"rot_index={rot_index}; cannot determine chirality sign"
            )
        if chirality_sign not in (1, -1, 1.0, -1.0):
            raise ValueError(
                f"chirality_index={chirality_index} requires chirality_sign "
                f"to be explicitly +1 or -1; got {chirality_sign!r}. This "
                "must be determined from the actual as-built array, not "
                "guessed or derived from the initial guess — see "
                "doc/SYMMETRY.md."
            )
        expected_chirality_sign = float(chirality_sign)
        chirality_norm = radius[rot_index] * radius[chirality_index]

    return _Problem(
        n_ant=n_ant,
        x0=x0,
        lower=lower,
        upper=upper,
        center=center,
        radius=np.asarray(radius, dtype=float),
        pairs=pairs,
        pair_values=pair_values,
        sqrt_radius_weight=np.sqrt(radius_weight),
        sqrt_rotation_weight=np.sqrt(rotation_weight),
        sqrt_chirality_weight=np.sqrt(chirality_weight),
        rot_index=rot_index,
        rot_degrees=rot_degrees,
        chirality_index=chirality_index,
        expected_chirality_sign=expected_chirality_sign,
        chirality_norm=chirality_norm,
    )


def _finalize_result(result, residuals_fn):
    """Normalise a least_squares result for interchange with callers that
    expect a scalar objective and an iteration count.

    ``scipy.optimize.least_squares`` stores the residual *vector* in
    ``result.fun`` and exposes no ``.nit``.  We overwrite ``.fun`` with the
    scalar sum-of-squares objective and add ``.nit`` (Jacobian evaluations
    as a proxy) so the result matches the historical contract of
    :func:`calibrate` / :func:`calibrate_irls`.
    """
    result.fun = float(np.sum(residuals_fn(result.x) ** 2))
    result.nit = getattr(result, "njev", getattr(result, "nfev", 0))
    return result


def calibrate(
    radius,
    m_ij,
    initial_guess,
    *,
    center=(0.0, 0.0),
    rot_index=None,
    rot_degrees=None,
    chirality_index=None,
    chirality_sign=None,
    chirality_weight=1e6,
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
        Only the upper triangle (i < j) is used so each pair is counted once.
    initial_guess : (N, 2) array
        Initial antenna positions in mm.
    center : (float, float)
        Reference point coordinates (mm).
    rot_index : int or None
        Antenna index used to constrain global rotation.  None disables.
    rot_degrees : float or None
        Target geographic angle (degrees) for *rot_index*.
    chirality_index : int or None
        Antenna index used to break the reflection (chirality) degeneracy.
        Requires *rot_index* and *chirality_sign*.
    chirality_sign : {1, -1, None}
        The required sign of the cross product ``p_ref × p_chirality``,
        enforced as a soft constraint.  This is **not** derived from
        *initial_guess* — it must be determined by the caller from the
        actual, as-built antenna array (e.g. by observing whether the
        chirality antenna lies to the left (+1) or right (-1) of the
        reference antenna's bearing line).  An as-built array can be
        mirrored relative to its nominal design file, so trusting the
        initial guess's sign is not safe; see ``doc/SYMMETRY.md``.
        Required (and validated) whenever *chirality_index* is given.
    chirality_weight : float
        Weight applied to the chirality penalty term.  Normalized cross
        product is dimensionless (≈sin(angle)), so a weight of 1e6
        effectively locks the chirality to *chirality_sign*.
    max_position_error : float
        Half-width of the search bounds around each initial coordinate (mm).
    radius_weight : float
        Weight applied to the radius residual term.
    rotation_weight : float
        Weight applied to the rotation constraint term.
    maxiter : int
        Maximum number of function evaluations for the optimiser.
    tol : float or None
        Convergence tolerance — sets ``ftol``, ``xtol`` and ``gtol`` of the
        trust-region (``trf``) least-squares solver.  ``None`` uses the
        scipy defaults.

    Returns
    -------
    result : OptimizeResult
        Result of :func:`scipy.optimize.least_squares` (method ``"trf"``).
        For convenience ``result.fun`` is overwritten with the scalar
        sum-of-squares objective and ``result.nit`` holds the Jacobian
        evaluation count.  ``result.x`` is the flat parameter vector
        ``(x0, y0, x1, y1, …)`` in mm.

    Raises
    ------
    ValueError
        If ``initial_guess`` or ``m_ij`` dimensions disagree with
        ``len(radius)``, if ``rot_index``/``chirality_index`` are outside
        ``[0, n_ant)``, or if ``chirality_index`` is given without a valid
        ``chirality_sign`` of ``+1`` or ``-1``.
    """
    prob = _build_problem(
        radius,
        m_ij,
        initial_guess,
        center=center,
        rot_index=rot_index,
        rot_degrees=rot_degrees,
        chirality_index=chirality_index,
        chirality_sign=chirality_sign,
        chirality_weight=chirality_weight,
        max_position_error=max_position_error,
        radius_weight=radius_weight,
        rotation_weight=rotation_weight,
    )

    if tol is not None:
        result = least_squares(
            prob.combined_residuals,
            prob.x0,
            bounds=prob.bounds,
            method="trf",
            max_nfev=maxiter,
            ftol=tol,
            xtol=tol,
            gtol=tol,
        )
    else:
        result = least_squares(
            prob.combined_residuals,
            prob.x0,
            bounds=prob.bounds,
            method="trf",
            max_nfev=maxiter,
        )
    return _finalize_result(result, prob.combined_residuals)


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
    chirality_index=None,
    chirality_sign=None,
    chirality_weight=1e6,
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

    prob = _build_problem(
        radius,
        m_ij,
        initial_guess,
        center=center,
        rot_index=rot_index,
        rot_degrees=rot_degrees,
        chirality_index=chirality_index,
        chirality_sign=chirality_sign,
        chirality_weight=chirality_weight,
        max_position_error=max_position_error,
        radius_weight=radius_weight,
        rotation_weight=rotation_weight,
    )

    n_meas = prob.n_radius + prob.n_ij
    weights = np.ones(prob.n_total)
    prev_x = prob.x0.copy()

    for _ in range(irls_max_iter):

        def weighted_residuals(x):
            return np.sqrt(weights) * prob.combined_residuals(x)

        res_ls = least_squares(
            weighted_residuals,
            prev_x,
            bounds=prob.bounds,
            method="trf",
            max_nfev=maxiter,
            xtol=1e-12,
            ftol=1e-12,
        )
        new_x = res_ls.x

        # Update weights from measurement residuals only (rotation and
        # chirality are priors, not noisy data, so they stay at weight 1).
        raw_res = prob.combined_residuals(new_x)
        weights[:n_meas] = weight_fn(raw_res[:n_meas])

        shift = np.max(np.abs(new_x - prev_x))
        prev_x = new_x

        if shift < irls_tol:
            break

    res_ls.x = prev_x
    return _finalize_result(res_ls, prob.combined_residuals)


def result_to_positions(result, n_ant):
    """Extract (n_ant, 2) positions array from an OptimizeResult."""
    return result.x.reshape((n_ant, 2))


def result_to_json(result, n_ant, title=None):
    """Convert an OptimizeResult to the API-compatible JSON dict.

    Positions are rounded to 3 decimal places (metres) with z=0.
    If *title* is provided, it is included as a ``"title"`` key.
    """
    pos_mm = result_to_positions(result, n_ant)
    out = np.zeros((n_ant, 3))
    out[:, :2] = np.round(pos_mm / 1000.0, 3)
    d = {"antenna_positions": out.tolist()}
    if title is not None:
        d["title"] = title
    return d


def compute_residuals(result, radius, m_ij, n_ant):
    """Return (radius_residuals, ij_residuals, ij_pairs) for diagnostics.

    Inter-antenna pairs are reported once each (upper triangle, i < j).
    """
    x = result.x

    r_res = np.array([dist((0, 0), _p(x, i)) - radius[i] for i in range(n_ant)])

    pairs = []
    ij_res = []
    for i in range(n_ant):
        for j in range(i + 1, n_ant):
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

    abs_ij = np.abs(ij_res)
    mad = np.median(np.abs(ij_res - np.median(ij_res)))
    std = np.std(ij_res)
    p50 = np.percentile(abs_ij, 50)
    p90 = np.percentile(abs_ij, 90)

    print("\nInter-antenna residual statistics:")
    print(f"  Median Absolute Deviation: {mad:.2f} mm")
    print(f"  Standard Deviation:       {std:.2f} mm")
    print(f"  50th percentile |res|:    {p50:.2f} mm")
    print(f"  90th percentile |res|:    {p90:.2f} mm")

    big = []
    for r, (i, j) in zip(ij_res, pairs):
        if np.abs(r) > p90:
            big.append((r, i, j))

    if big:
        print("\nLargest residuals (>90th pct):")
        for r, i, j in sorted(big):
            print(f"  res[{i},{j}] = {r:+.1f}")
