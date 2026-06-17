"""CLI entry point for tart-position-cal antenna position calibration."""

import argparse
import json
import sys

from tart_position_cal.calibrate import (
    calibrate,
    calibrate_irls,
    compute_residuals,
    fetch_initial_positions,
    load_initial_positions,
    load_measurements,
    plot_differences,
    plot_positions,
    plot_residual_histogram,
    print_residual_report,
    result_to_json,
    result_to_positions,
)


def main():
    parser = argparse.ArgumentParser(
        description="Calibrate TART antenna positions from survey measurements."
    )
    parser.add_argument(
        "--measurements",
        required=True,
        help="Path to .ods file with antenna distance measurements.",
    )
    parser.add_argument(
        "--positions",
        help="Path to JSON file with initial antenna positions (metres).",
    )
    parser.add_argument(
        "--api",
        help="TART API base URL to fetch initial positions from (e.g. https://api.elec.ac.nz/tart/na-unam).",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path for output calibrated positions JSON file.",
    )
    parser.add_argument(
        "--rot-index",
        type=int,
        default=0,
        help="Antenna index used for global rotation constraint.",
    )
    parser.add_argument(
        "--rot-degrees",
        type=float,
        default=0.0,
        help="Target geographic angle in degrees for the constrained antenna.",
    )
    parser.add_argument(
        "--max-iter",
        type=int,
        default=500,
        help="Maximum optimizer iterations (default: 500).",
    )
    parser.add_argument(
        "--max-position-error",
        type=float,
        default=4200.0,
        help="Half-width of search bounds around initial positions in mm (default: 4200).",
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Skip generation of diagnostic plots.",
    )
    parser.add_argument(
        "--plot-dir",
        default=".",
        help="Directory to save diagnostic plots (default: current directory).",
    )
    parser.add_argument(
        "--irls",
        action="store_true",
        help="Use Iteratively Reweighted Least Squares for robust calibration.",
    )
    parser.add_argument(
        "--irls-weight-fn",
        choices=["tukey", "huber"],
        default="tukey",
        help="Weight function for IRLS: 'tukey' (bisquare) or 'huber' (default: tukey).",
    )
    parser.add_argument(
        "--irls-max-iter",
        type=int,
        default=10,
        help="Maximum IRLS outer iterations (default: 10).",
    )

    args = parser.parse_args()

    # --- Load measurements ---
    print(f"Loading measurements from {args.measurements}")
    radius, m_ij, n_ant = load_measurements(args.measurements)
    print(f"  Found {n_ant} antennas")

    # --- Load initial positions ---
    if args.positions:
        print(f"Loading initial positions from {args.positions}")
        initial_guess = load_initial_positions(args.positions)
    elif args.api:
        print(f"Fetching initial positions from {args.api}")
        initial_guess = fetch_initial_positions(args.api)
    else:
        parser.error("One of --positions or --api is required.")
        sys.exit(1)

    if initial_guess.shape[0] != n_ant:
        print(
            f"WARNING: measurements have {n_ant} antennas but "
            f"positions have {initial_guess.shape[0]}"
        )

    # --- Calibrate ---
    if args.irls:
        print(
            f"Running IRLS calibration (rot_index={args.rot_index}, "
            f"rot_degrees={args.rot_degrees}, "
            f"weight_fn={args.irls_weight_fn})"
        )
        result = calibrate_irls(
            radius,
            m_ij,
            initial_guess,
            rot_index=args.rot_index,
            rot_degrees=args.rot_degrees,
            max_position_error=args.max_position_error,
            maxiter=args.max_iter,
            irls_max_iter=args.irls_max_iter,
            weight_function=args.irls_weight_fn,
        )
        print(f"  Final objective: {result.fun:.4f}")
        print(f"  Function evaluations: {result.nfev}")
    else:
        print(
            f"Running calibration (rot_index={args.rot_index}, "
            f"rot_degrees={args.rot_degrees})"
        )
        result = calibrate(
            radius,
            m_ij,
            initial_guess,
            rot_index=args.rot_index,
            rot_degrees=args.rot_degrees,
            max_position_error=args.max_position_error,
            maxiter=args.max_iter,
        )

        print(f"  Optimisation finished: {result.message}")
        print(f"  Final objective: {result.fun:.4f}")
        print(f"  Iterations: {result.nit}")
        print(f"  Function evaluations: {result.nfev}")

    # --- Save output ---
    json_out = result_to_json(result, n_ant)
    with open(args.output, "w") as f:
        json.dump(json_out, f)
    print(f"Calibrated positions written to {args.output}")

    # --- Diagnostics ---
    print_residual_report(result, radius, m_ij, n_ant)

    # --- Plots ---
    if not args.no_plots:
        import os

        final_pos = result_to_positions(result, n_ant)
        plot_positions(
            initial_guess,
            final_pos,
            os.path.join(args.plot_dir, "final_positions.png"),
        )
        print(f"Saved final_positions.png to {args.plot_dir}")

        plot_differences(
            initial_guess,
            final_pos,
            os.path.join(args.plot_dir, "differences.png"),
        )
        print(f"Saved differences.png to {args.plot_dir}")

        _, ij_res, _ = compute_residuals(result, radius, m_ij, n_ant)
        plot_residual_histogram(
            ij_res,
            os.path.join(args.plot_dir, "residual_histogram.png"),
        )
        print(f"Saved residual_histogram.png to {args.plot_dir}")
