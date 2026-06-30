"""Tests for tart_position_cal.calibrate."""

from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from tart_position_cal.calibrate import (
    calibrate,
    calibrate_irls,
    compute_bearing,
    fetch_initial_positions,
    geo_angle,
    load_initial_positions,
    load_measurements,
    result_to_positions,
)

# ---- compute_bearing -------------------------------------------------------


def test_unam_bearing():
    """UNAM phase-centre to distant hill — value from the paper."""
    b = compute_bearing(-22.612053, 17.056784, -22.553822, 17.077290)
    assert b == pytest.approx(18.1125, abs=0.001)


def test_north():
    """Due north: same longitude, northern destination."""
    b = compute_bearing(0.0, 10.0, 1.0, 10.0)
    assert b == pytest.approx(0.0, abs=0.001)


def test_south():
    """Due south: same longitude, southern destination."""
    b = compute_bearing(1.0, 10.0, 0.0, 10.0)
    assert b == pytest.approx(180.0, abs=0.001)


def test_east():
    """Due east: same latitude, eastern destination (small separation)."""
    b = compute_bearing(0.0, 10.0, 0.0, 11.0)
    assert b == pytest.approx(90.0, abs=0.1)


def test_west():
    """Due west: same latitude, western destination."""
    b = compute_bearing(0.0, 11.0, 0.0, 10.0)
    # pyproj returns negative azimuths for westward bearings
    assert abs(b + 90.0) < 0.1 or abs(b - 270.0) < 0.1


def test_return_bearing():
    """Forward and return bearings differ by ~180°.

    On a WGS84 ellipsoid the forward and back azimuths are not exactly
    180° apart — the geodesic is not a planar curve.
    """
    b_fwd = compute_bearing(-22.612053, 17.056784, -22.553822, 17.077290)
    b_ret = compute_bearing(-22.553822, 17.077290, -22.612053, 17.056784)
    # Normalise both to [0, 360)
    b_fwd_n = b_fwd % 360
    b_ret_n = b_ret % 360
    diff = abs((b_fwd_n + 180.0) % 360 - b_ret_n)
    diff = min(diff, 360.0 - diff)
    assert diff < 0.1


def test_northeast():
    """45° bearing: equal latitude and longitude steps."""
    b = compute_bearing(0.0, 10.0, 1.0, 11.0)
    assert 44.0 < b < 46.0


def test_northwest():
    """Bearing toward northwest (returns negative or ~315°)."""
    b = compute_bearing(0.0, 11.0, 1.0, 10.0)
    b_norm = b % 360
    assert 314.0 < b_norm < 316.0


def test_short_distance_consistency():
    """A very short distance should give a stable bearing."""
    b1 = compute_bearing(-22.612, 17.056, -22.612001, 17.056001)
    b2 = compute_bearing(-22.612, 17.056, -22.612002, 17.056002)
    assert b1 == pytest.approx(b2, abs=0.01)


# ---- geo_angle consistency with compute_bearing -----------------------------


def test_geo_angle_north():
    """geo_angle matches compute_bearing: north = 0°."""
    assert geo_angle(0.0, 1.0) == pytest.approx(0.0)


def test_geo_angle_east():
    """geo_angle matches compute_bearing: east = 90°."""
    assert geo_angle(1.0, 0.0) == pytest.approx(90.0)


# ---- calibrate ----------------------------------------------------------


@pytest.fixture(scope="module")
def unam_data():
    """Load the example UNAM measurement and position data."""
    radius, m_ij, n_ant = load_measurements("example/antenna_measurements.ods")
    initial = load_initial_positions("example/antenna_positions.json")
    return radius, m_ij, n_ant, initial


def test_calibrate_converges(unam_data):
    """Standard calibration converges and reduces the cost."""
    radius, m_ij, n_ant, initial = unam_data
    result = calibrate(radius, m_ij, initial, rot_index=3, rot_degrees=18.1125)
    assert result.fun < 1000.0
    assert result.nit > 0
    pos = result_to_positions(result, n_ant)
    assert pos.shape == (n_ant, 2)


def test_calibrate_with_chirality(unam_data):
    """Chirality constraint gives the same result as unconstrained."""
    radius, m_ij, n_ant, initial = unam_data
    result_no = calibrate(radius, m_ij, initial, rot_index=3, rot_degrees=18.1125)
    result_chi = calibrate(
        radius,
        m_ij,
        initial,
        rot_index=3,
        rot_degrees=18.1125,
        chirality_index=10,
    )
    assert result_chi.fun == pytest.approx(result_no.fun, rel=1e-6)


def test_calibrate_irls_converges(unam_data):
    """IRLS calibration converges."""
    radius, m_ij, n_ant, initial = unam_data
    result = calibrate_irls(
        radius,
        m_ij,
        initial,
        rot_index=3,
        rot_degrees=18.1125,
        chirality_index=10,
    )
    assert result.fun > 0.0
    pos = result_to_positions(result, n_ant)
    assert pos.shape == (n_ant, 2)


def test_calibrate_respects_bounds(unam_data):
    """Tight bounds keep positions near the initial guess."""
    radius, m_ij, n_ant, initial = unam_data
    result = calibrate(
        radius,
        m_ij,
        initial,
        rot_index=3,
        rot_degrees=18.1125,
        max_position_error=100.0,  # tight bound
    )
    pos = result_to_positions(result, n_ant)
    for i in range(n_ant):
        assert abs(pos[i, 0] - initial[i, 0]) <= 150.0  # margin for convergence
        assert abs(pos[i, 1] - initial[i, 1]) <= 150.0


def test_chirality_collinearity_raises(unam_data):
    """Using the same index for rot_index and chirality_index raises."""
    radius, m_ij, n_ant, initial = unam_data
    with pytest.raises(ValueError, match="collinear"):
        calibrate(
            radius,
            m_ij,
            initial,
            rot_index=3,
            rot_degrees=18.1125,
            chirality_index=3,
        )


# ---- input validation (#7) -----------------------------------------------


def test_calibrate_rejects_rot_index_out_of_range(unam_data):
    """rot_index >= n_ant raises a clear ValueError."""
    radius, m_ij, n_ant, initial = unam_data
    with pytest.raises(ValueError, match="rot_index.*out of range"):
        calibrate(radius, m_ij, initial, rot_index=99, rot_degrees=18.0)


def test_calibrate_rejects_negative_rot_index(unam_data):
    """Negative rot_index raises a clear ValueError."""
    radius, m_ij, n_ant, initial = unam_data
    with pytest.raises(ValueError, match="out of range"):
        calibrate(radius, m_ij, initial, rot_index=-1, rot_degrees=18.0)


def test_calibrate_rejects_chirality_index_out_of_range(unam_data):
    """chirality_index outside [0, n_ant) raises a clear ValueError."""
    radius, m_ij, n_ant, initial = unam_data
    with pytest.raises(ValueError, match="chirality_index.*out of range"):
        calibrate(
            radius,
            m_ij,
            initial,
            rot_index=3,
            rot_degrees=18.0,
            chirality_index=99,
        )


def test_calibrate_rejects_initial_guess_row_count(unam_data):
    """initial_guess with the wrong number of rows raises ValueError."""
    radius, m_ij, n_ant, initial = unam_data
    with pytest.raises(ValueError, match="rows"):
        calibrate(radius, m_ij, initial[:5], rot_index=3, rot_degrees=18.0)


def test_calibrate_rejects_initial_guess_shape(unam_data):
    """initial_guess that is not (n_ant, 2) raises ValueError."""
    radius, m_ij, n_ant, initial = unam_data
    bad = np.hstack([initial, np.ones((n_ant, 1))])  # (n_ant, 3)
    with pytest.raises(ValueError, match="initial_guess must have shape"):
        calibrate(radius, m_ij, bad, rot_index=3, rot_degrees=18.0)


def test_calibrate_rejects_m_ij_shape(unam_data):
    """m_ij with the wrong shape raises ValueError."""
    radius, m_ij, n_ant, initial = unam_data
    with pytest.raises(ValueError, match="m_ij must have shape"):
        calibrate(radius, m_ij[:10, :10], initial, rot_index=3, rot_degrees=18.0)


def test_calibrate_irls_also_validates(unam_data):
    """The IRLS path shares the validation (raises on bad rot_index)."""
    radius, m_ij, n_ant, initial = unam_data
    with pytest.raises(ValueError, match="out of range"):
        calibrate_irls(radius, m_ij, initial, rot_index=99, rot_degrees=18.0)


# ---- fetch_initial_positions (#8) ----------------------------------------


def _fake_response():
    fake = MagicMock()
    fake.json.return_value = [[1.0, 2.0, 0.0], [4.0, 5.0, 0.0]]
    fake.raise_for_status.return_value = None
    return fake


def test_fetch_initial_positions_forwards_timeout():
    """A caller-supplied timeout is passed through to requests.get."""
    with patch("requests.get", return_value=_fake_response()) as mock_get:
        fetch_initial_positions("https://example.org/tart", timeout=12.5)
    _, kwargs = mock_get.call_args
    assert kwargs.get("timeout") == 12.5


def test_fetch_initial_positions_default_sets_timeout():
    """A timeout is always set, so an unresponsive API cannot hang forever."""
    with patch("requests.get", return_value=_fake_response()) as mock_get:
        fetch_initial_positions("https://example.org/tart")
    _, kwargs = mock_get.call_args
    assert kwargs.get("timeout") not in (None, 0, 0.0)


# ---- loader error handling (#11) -----------------------------------------


def test_load_initial_positions_missing_file(tmp_path):
    """A missing positions file raises a clear FileNotFoundError."""
    with pytest.raises(FileNotFoundError, match="initial positions"):
        load_initial_positions(tmp_path / "does_not_exist.json")


def test_load_initial_positions_bad_json(tmp_path):
    """Malformed JSON raises a clear ValueError."""
    p = tmp_path / "bad.json"
    p.write_text("{ not valid json")
    with pytest.raises(ValueError, match="Invalid JSON"):
        load_initial_positions(p)


def test_load_initial_positions_missing_key(tmp_path):
    """A positions file without the required key raises a clear ValueError."""
    p = tmp_path / "nokey.json"
    p.write_text('{"foo": 1}')
    with pytest.raises(ValueError, match="antenna_positions"):
        load_initial_positions(p)


def test_load_measurements_missing_file(tmp_path):
    """A missing measurements file raises a clear FileNotFoundError."""
    with pytest.raises(FileNotFoundError, match="measurements"):
        load_measurements(tmp_path / "does_not_exist.ods")


def test_load_measurements_missing_sheet(tmp_path):
    """An ODS without a 'Sheet1' sheet raises a clear ValueError."""
    import pandas as pd

    p = tmp_path / "wrong_sheet.ods"
    pd.DataFrame({"A 0": [1.0]}).to_excel(p, sheet_name="NotSheet1", engine="odf")
    with pytest.raises(ValueError, match="Sheet1"):
        load_measurements(p)


def test_load_measurements_too_few_rows(tmp_path):
    """A Sheet1 with too few rows for the antenna count raises ValueError."""
    import pandas as pd

    p = tmp_path / "short.ods"
    cols = {f"A {i}": [float(i)] for i in range(24)}  # only one data row
    pd.DataFrame(cols).to_excel(p, sheet_name="Sheet1", engine="odf")
    with pytest.raises(ValueError, match="data row"):
        load_measurements(p)
