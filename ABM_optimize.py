import argparse
import csv
import json
import os
import shutil
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from ABM import (
    create_sample_file,
    run_ABM,
    extract_output_metrics,
    metric_label,
    CONFIG_TEMPLATE,
    SAMPLE_CONFIG_PATH,
    OUTPUT_BIOMARKERS,
    OUTPUT_DIR,
    PARAMETER_FILE,
    CONDITIONS,
    TICKS_PER_DAY,
    OUTPUT_METRICS,
    FITTED_METRICS,
    FITTED_METRIC_KEYS,
    METRICS_NUMTICKS,
)

# ==========================
# CONSTANTS
# TODO: Move these into a config file.
SNAPSHOT_INTERVAL = 5

CONDITION_TO_GROUP = {
    "high": "config_scaffold_High",
    "low": "config_scaffold_Low",
}

# Objective function's error terms
# These come from OUTPUT_METRICS in ABM.py, specifically those carrying
# an "exp_column" (i.e. there is exp data for that metric that we can fit
# against).
TARGETS = FITTED_METRICS
 
# Simulate far enough to cover the latest day any metric needs
NUMTICKS = METRICS_NUMTICKS
 
SENSITIVITY_DIR = OUTPUT_DIR / "SensitivityAnalysis"
# ==========================

# ==========================
# Command line arguments

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run ABM Nelder-Mead optimization.")
 
    parser.add_argument(
        "--mode", choices=["joint", "separate"], default="joint",
        help=(
            "How to handle the scaffold conditions. 'joint' (default) runs one "
            "optimization whose objective is the error summed over every "
            "condition, giving a single parameter set that fits all of them. "
            "'separate' runs an independent optimization per condition, giving "
            "one parameter set and one error per condition"
        )
    )
    parser.add_argument(
        "--conditions", nargs="+", default=CONDITIONS, choices=CONDITIONS,
        help="Which scaffold conditions to fit (default: all)"
    )
    parser.add_argument(
        "--iters", type=int, default=3,
        help="ABM runs per condition per function evaluation, averaged"
    )
    parser.add_argument(
        "--maxiter", type=int, default=50,
        help="Maximum Nelder-Mead iterations"
    )
    parser.add_argument(
        "--tol", type=float, default=1e-4,
        help="Nelder-Mead convergence tolerance"
    )
    parser.add_argument(
        "--snapshots", action="store_true",
        help="Save periodic snapshots of the biomarker values"
    )
    parser.add_argument(
        "--out", type=Path, default=OUTPUT_DIR / "optimization_results.json",
        help="Where to write the optimization results"
    )
    return parser.parse_args(argv)
# ==========================

def small_scaffold_adjustment_cells(value: float) -> float:
    """
    Adjusts the value for the small scaffold experimental data.
    The adjustment is based on the formula:
    ((0.6)^3 / 1000 / 0.3 ) * value
    """
    return ((0.3 ** 3) / 1000 / 0.3) * value

def small_scaffold_adjustment_aggrecan(value: float) -> float:
    """
    Adjusts the value for the small scaffold experimental data.
    The adjustment is based on the formula:
    ((0.6)^3 / 300 ) * value * 10^6
    """
    return ((0.5 * 0.4 * 0.3) / 60) * value * 10 ** 6 # value in ug

def extract_small_scaffold_experimental(file_path: Path) -> pd.DataFrame:
    """
    Extracts experimental data from CSV file corresponding to the small scaffold.
    """
    df = pd.read_csv(file_path)
    df.columns = [c.strip() for c in df.columns]

    # Calculate average for each config (low, high)
    df["cell_viability_percent"] = df["cell_viability_percent"].astype(float)
    df["sGAG_total_ug"] = df["sGAG_total_ug"].astype(float)
    df["percent_diff"] = df["percent_diff"].astype(float)

    # Calculate the average values for each group and time point
    averages = df.groupby(["group", "time_hour"]).mean().reset_index()
    
    # Cell-related adjustments
    averages["small_scaffold_cell_viability"] = averages["cell_viability_percent"]
    averages["small_scaffold_percent_diff"] = averages["percent_diff"]
    
    # Aggrecan adjustment
    averages["small_scaffold_aggrecan_ug"] = averages["sGAG_total_ug"].apply(
        lambda x: small_scaffold_adjustment_aggrecan(x)
    )
    
    return averages

def error(expected, actual):
    """
    Calculate SSE between expected and actual values.
    """
    return (expected - actual) ** 2

def construct_simplex(bounds: np.ndarray, selected_params: list):
    """
    Construct a simplex for n selected parameters
    by perturbing the default values (beginning with the minimums)
    changing one to a maximum at a time.

    The simplex consists of n+1 points in n-dimensional space.
    Each point is represented by a list of n values.

    Example:
    init_simplex = [
        [2,0,0,2.5,1],
        [200, 0, 0,2.5,1],
        [2, 50, 0,2.5,1],
        [2, 0, 9630,2.5,1],
        [2,0,0,100,1],
        [200,50,9639,100,100]
    ]
    """
    # create list with just the minimums to start
    # and store the list in init_simplex
    init_simplex = [[bounds[selected_param][0] for selected_param in selected_params]]

    # selectively choose one parameter at a time
    # to set to the max value
    for i in range(len(selected_params) - 1):
        # create a copy of the simplex
        simplex = init_simplex[0].copy()
        # set the selected parameter to the maximum
        simplex[i] = bounds[selected_params[i]][1]
        # append the new point to the simplex
        init_simplex.append(simplex)
    
    # in the last row, set all parameters to the maximum
    # and append to the simplex
    simplex = init_simplex[0].copy()
    for i in range(len(selected_params)):
        simplex[i] = bounds[selected_params[i]][1]
    init_simplex.append(simplex)

    print(f"""
        Initial simplex:
        {[[float(value) for value in row] for row in init_simplex]}
        """)
    
    return init_simplex

def formatted_string(Nfeval, x, Y) -> str:
    """
    Create dynamic string in this format representing the parameters:
    Nfeval x1 ... xn+1 SSE
    With the following format:
    {0:4d} {1: 3.6f} ... {n+1: 3.6f} {error: 3.6f}
    Where n is the number of parameters
    """
    # iteration number
    format_str = f"{Nfeval:4d} "
    # add parameter with index
    for i, param in enumerate(x, start=1):
        format_str += f"{i}: {param:3.6f} "
    # add sum of Y at end
    format_str += f"{np.sum(Y):3.6f}"
    return format_str

def extract_varying_param_indices(df: pd.DataFrame) -> list:
    """
    Extract the indices of every parameter that should be optimized --
    i.e. every row in parameters.xlsx whose "Vary?" column is
    blank/NaN, matching the same convention used in the upstream
    sampling/RF pipeline (blank = vary, any non-blank value e.g. "N" =
    don't vary, use the default).
    """
    vary_mask = df["Vary?"].apply(lambda v: not (pd.notna(v) and str(v).strip() != ""))
    varying_rows = df[vary_mask]
 
    if varying_rows.empty:
        raise ValueError('No parameters found with a blank "Vary?" column -- nothing to optimize.')
 
    param_idxs = [int(p) - 1 for p in varying_rows["Parameter Number"].tolist()]

    if any(i < 0 or i >= len(df) for i in param_idxs):
        raise ValueError(
            '"Parameter Number" values in parameters.xlsx are out of range; '
            "they must be 1-indexed and contiguous with the parameter rows."
        )

    print(f"Optimizing {len(param_idxs)} parameters (Vary? blank), rows: "
          f"{[i + 1 for i in param_idxs]}")

    return param_idxs

def extract_row_as_dict(csv_path: str, row_index: int) -> dict[str, Any]:
    """
    Read the OutputBiomarkers.csv file which is the output of a single run
    of the ABM simulation.
    Extract all biomarker values from the specified row.
    """
    with open(csv_path, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for i, row in enumerate(reader):
            if i == row_index:
                return {k: float(v) for k, v in row.items()}
    return {}

def save_snapshot(nfeval: int, condition: str):
    """
    Save snapshot for this condition to a CSV file.
    """
    rows = []
    for day in sorted({m["day"] for m in OUTPUT_METRICS}):
        snapshot = extract_row_as_dict(str(OUTPUT_BIOMARKERS), TICKS_PER_DAY * day)
        if not snapshot:
            continue
        snapshot["day"] = day
        snapshot["nfeval"] = nfeval
        snapshot["condition"] = condition
        rows.append(snapshot)
 
    if not rows:
        return
 
    fieldnames = list(rows[0].keys())
    snapshot_path = OUTPUT_DIR / "snapshots" / "day_snapshots.csv"
    with open(snapshot_path, 'a', newline='') as snapshots_file:
        writer = csv.DictWriter(snapshots_file, fieldnames=fieldnames)
        if snapshots_file.tell() == 0:
            writer.writeheader()
        writer.writerows(rows)

def run_with_scaffold(
    condition: str,
    experimental_df: pd.DataFrame,
    param_values: list,
    param_names: list,
    nfeval: int,
    num_iters: int = 3,
) -> dict[str, float]:
    """
    Run the ABM with the specified scaffold condition ("high" or "low")
    num_iters times, and return the mean squared error per fitted target,
    keyed by the target's label.
 
    Every metric in OUTPUT_METRICS is read and printed each run, but only
    those with an "exp_column" (i.e. TARGETS) produce error terms
    """
    group_name = CONDITION_TO_GROUP[condition]
    errors_per_target: dict[str, list[float]] = {
        metric_label(t): [] for t in TARGETS
    }
 
    for iteration in range(num_iters):
        print(f"Running iteration {iteration + 1}/{num_iters} for condition '{condition}'")
 
        with open(SENSITIVITY_DIR / "stdout.txt", 'a') as stdout_file, \
             open(SENSITIVITY_DIR / "stderr.txt", 'a') as stderr_file:
            banner = (
                "\n\n******************************\n"
                f"*** MODEL EXECUTION #{nfeval} ({condition}) ***\n"
                "******************************\n"
            )
            stdout_file.write(banner)
            stderr_file.write(banner)
 
        create_sample_file(
            param_values, CONFIG_TEMPLATE, SAMPLE_CONFIG_PATH,
            parameter_names=param_names, condition=condition,
        )
        run_ABM(SAMPLE_CONFIG_PATH, numticks=NUMTICKS)
 
        # Read every metric in OUTPUT_METRICS in one pass, keyed by label
        simulated = extract_output_metrics(OUTPUT_BIOMARKERS)
        print(f"  {condition}: " + ", ".join(
            f"{label}={value:.4g}" for label, value in simulated.items()
        ))
 
        # Error terms come from the fitted subset only
        for target in TARGETS:
            label = metric_label(target)
            expected = float(
                experimental_df.loc[(group_name, target["day"] * 24), target["exp_column"]]
            )
            errors_per_target[label].append(error(expected, simulated[label]))
 
    # Mean error over all runs
    return {label: float(np.mean(errs)) for label, errs in errors_per_target.items()}

def make_objective(
    conditions: list[str],
    params: list[int],
    default_values: np.ndarray,
    param_names: list[str],
    experimental_df: pd.DataFrame,
    args: argparse.Namespace,
) -> Callable[[np.ndarray], float]:
    """
    Build the Nelder-Mead objective function.
 
    The returned function evaluates the ABM once per condition in
    `conditions` and returns the total SSE summed over all of them. Pass a
    single condition for a per-condition fit, or every condition for a
    joint fit -- the maths is the same, only the length of `conditions`
    differs.
    """
    state = {"nfeval": 1, "history": []}
 
    def objective(x: np.ndarray) -> float:
        nfeval = state["nfeval"]
 
        # Start from the spreadsheet defaults and overwrite only the
        # parameters being optimized with the simplex's current values
        sam = np.array(default_values, dtype=float)
        for idx, param in enumerate(params):
            sam[param] = x[idx]
        param_values = sam.tolist()
 
        errors_by_condition: dict[str, dict[str, float]] = {}
 
        for condition in conditions:
            errors_by_condition[condition] = run_with_scaffold(
                condition,
                experimental_df,
                param_values=param_values,
                param_names=param_names,
                nfeval=nfeval,
                num_iters=args.iters,
            )
 
            if args.snapshots:
                save_snapshot(nfeval, condition)
                if nfeval % SNAPSHOT_INTERVAL == 0:
                    shutil.copy(
                        OUTPUT_BIOMARKERS,
                        OUTPUT_DIR / "snapshots" / "biomarker_csvs" /
                        f"snapshots_{nfeval}_{condition}.csv",
                    )
 
        # Flatten to the objective vector Y: one entry per target per
        # condition, in (condition, target) order
        Y = np.array([
            errors_by_condition[condition][metric_label(target)]
            for condition in conditions
            for target in TARGETS
        ])
 
        total = float(np.sum(Y))
 
        print(formatted_string(nfeval, x, Y))
        print({c: errors_by_condition[c] for c in conditions})
        print(f"  total SSE over {len(conditions)} condition(s) "
              f"({', '.join(conditions)}): {total:.6f}")
 
        state["history"].append({
            "nfeval": nfeval,
            "x": [float(v) for v in x],
            "errors_by_condition": errors_by_condition,
            "total_sse": total,
        })
        state["nfeval"] = nfeval + 1
 
        return total  # SSE
 
    objective.state = state
    return objective

def run_optimization(
    conditions: list[str],
    params: list[int],
    bounds: np.ndarray,
    default_values: np.ndarray,
    param_names: list[str],
    experimental_df: pd.DataFrame,
    args: argparse.Namespace,
) -> dict:
    """
    Run one Nelder-Mead optimization fitting `conditions` together, and
    return a summary dict of the result
    """
    label = "+".join(conditions)
    print(f"\n{'=' * 70}\nOptimizing conditions: {label}\n{'=' * 70}")
 
    objective = make_objective(
        conditions, params, default_values, param_names, experimental_df, args
    )
 
    initial_simplex = construct_simplex(bounds, params)
    selected_bounds = bounds[np.array(params)]
    selected_defaults = default_values[np.array(params)]
 
    result = minimize(
        objective,
        selected_defaults,
        method='nelder-mead',
        tol=args.tol,
        bounds=selected_bounds,
        options={
            'maxiter': args.maxiter,
            'disp': True,
            'initial_simplex': initial_simplex,
            'return_all': True,
        }
    )
 
    history = objective.state["history"]
    initial_sse = history[0]["total_sse"] if history else None
    final_sse = float(result.fun)
 
    print(f"\nResults for {label}:")
    print(f"  optimal parameters: {result.x}")
    print(f"  minimum SSE: {final_sse}")
    print(f"  stopping criterion: {result.message}")
    if initial_sse:
        print(f"  SSE decrease: {initial_sse - final_sse:.6f} "
              f"({100 * (initial_sse - final_sse) / initial_sse:.2f}%)")
 
    # Full parameter vector with the optimized values written back in
    optimized_full = np.array(default_values, dtype=float)
    for idx, param in enumerate(params):
        optimized_full[param] = result.x[idx]
 
    return {
        "conditions": conditions,
        "optimized_parameters": {
            param_names[param]: float(result.x[idx])
            for idx, param in enumerate(params)
        },
        "optimized_parameter_vector": {
            name: float(value) for name, value in zip(param_names, optimized_full)
        },
        "initial_sse": initial_sse,
        "final_sse": final_sse,
        "sse_decrease": (initial_sse - final_sse) if initial_sse else None,
        "sse_percent_decrease": (
            100 * (initial_sse - final_sse) / initial_sse if initial_sse else None
        ),
        "final_errors_by_condition": history[-1]["errors_by_condition"] if history else {},
        "num_function_evaluations": len(history),
        "success": bool(result.success),
        "message": str(result.message),
        "history": history,
    }

if __name__ == "__main__":
    args = parse_args()
 
    # Read parameter definitions
    df = pd.read_excel(PARAMETER_FILE)
    numpar = len(df)
    param_names = df["Parameter Name"].tolist()
 
    # Parameters to optimize: every row with a blank "Vary?" column
    params = extract_varying_param_indices(df)
 
    # Read bounds and default (mean) values
    bounds = df[["Lower", "Upper"]].to_numpy()
    default_values = np.reshape(df[["Mean"]].to_numpy(), numpar)
 
    print("Selected parameters:")
    for i in params:
        print(f"  {i + 1}: {param_names[i]} "
              f"(default {default_values[i]}, bounds {tuple(bounds[i])})")
 
    # Prepare output directories / logs
    SENSITIVITY_DIR.mkdir(parents=True, exist_ok=True)
    open(SENSITIVITY_DIR / "stdout.txt", 'w').close()
    open(SENSITIVITY_DIR / "stderr.txt", 'w').close()
 
    if args.snapshots:
        (OUTPUT_DIR / "snapshots").mkdir(parents=True, exist_ok=True)
        (OUTPUT_DIR / "snapshots" / "biomarker_csvs").mkdir(parents=True, exist_ok=True)
        open(OUTPUT_DIR / "snapshots" / "day_snapshots.csv", 'w').close()
 
    # Generate experimental values
    experimental_df = extract_small_scaffold_experimental(Path("experimental_config.csv"))
    experimental_df_indexed = experimental_df.set_index(['group', 'time_hour'])
 
    print("Experimental data:")
    print(experimental_df_indexed)
 
    # # # Run optimization # # #
    # joint:    one fit, objective = error summed over every condition
    # separate: an independent fit (and parameter set) per condition
    if args.mode == "joint":
        condition_groups = [list(args.conditions)]
    else:
        condition_groups = [[condition] for condition in args.conditions]
 
    runs = [
        run_optimization(
            group, params, bounds, default_values,
            param_names, experimental_df_indexed, args,
        )
        for group in condition_groups
    ]
 
    results = {
        "mode": args.mode,
        "conditions": list(args.conditions),
        "settings": {
            "iters_per_evaluation": args.iters,
            "maxiter": args.maxiter,
            "tol": args.tol,
            "numticks": NUMTICKS,
            "fitted_targets": FITTED_METRIC_KEYS,
            "tracked_metrics": [metric_label(m) for m in OUTPUT_METRICS],
        },
        "runs": runs,
    }
 
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(results, f, indent=2)
 
    print(f"\n{'=' * 70}")
    print(f"Mode: {args.mode}")
    for run in runs:
        print(f"  {'+'.join(run['conditions'])}: final SSE {run['final_sse']:.6f}"
              + (f" ({run['sse_percent_decrease']:.2f}% decrease)"
                 if run["sse_percent_decrease"] is not None else ""))
    if args.mode == "separate" and len(runs) > 1:
        print(f"  sum of per-condition SSEs: {sum(r['final_sse'] for r in runs):.6f}")
    print(f"Wrote {args.out}")