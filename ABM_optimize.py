import argparse
import csv
import time
import hashlib
import json
import os
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from ABM import (
    create_sample_file,
    prepare_workdir,
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
    WORLD_SIZE_MM,
    WORLD_VOLUME_MM3,
    GEL_VOLUME_MM3,
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
        "--iters", type=int, default=2,
        help="ABM runs per condition per function evaluation, averaged"
    )
    parser.add_argument(
        "--maxiter", type=int, default=50,
        help="Maximum Nelder-Mead iterations (not ABM iterations)"
    )
    parser.add_argument(
        "--maxfev", type=int, default=None,
        help=(
            "Hard cap on objective evaluations -- the reliable way to bound "
            "a run's cost. The optimizer stops cleanly when it's reached and "
            "still writes results. Good for fitting a run inside a walltime."
        )
    )
    parser.add_argument(
        "--numticks", type=int, default=None,
        help=(
            "TESTING ONLY: override the simulated tick count. The default "
            "covers the latest day in OUTPUT_METRICS; a smaller value makes "
            "executions faster but will fail when a target's day isn't "
            "reached. Use with --max-day-override for a quick plumbing check."
        )
    )
    parser.add_argument(
        "--tol", type=float, default=1e-3,
        help="Nelder-Mead convergence tolerance"
    )
    parser.add_argument(
        "--parallel", type=int, default=None,
        help=(
            "How many ABM executions to run concurrently. Within one "
            "Nelder-Mead evaluation the (condition x --iters) executions are "
            "independent, so this is where the speedup is; Nelder-Mead itself "
            "is sequential and cannot be batched. Default: all of them at "
            "once (conditions x iters). Set 1 for the old serial behaviour."
            )
    )
    parser.add_argument(
        "--devices", type=int, default=None,
        help=(
            "Number of GPUs to spread concurrent executions across, "
            "round-robin via CUDA_VISIBLE_DEVICES. Default: read from "
            "SLURM_GPUS_ON_NODE, else 1."
        )
    )
    parser.add_argument(
        "--cache", type=Path, default=OUTPUT_DIR / "objective_cache.json",
        help=(
            "Checkpoint file of evaluated parameter vectors -> errors. Lets a "
            "timed-out run resume: Nelder-Mead is deterministic given the "
            "objective values, so replaying with cached results fast-forwards "
            "to where the previous job stopped."
        )
    )
    parser.add_argument(
        "--no-cache", action="store_true",
        help="Ignore and don't write the cache; every evaluation re-runs the ABM."
    )
    parser.add_argument(
        "--workdir-root", type=Path, default=None,
        help=(
            "Where per-execution working directories go. Default: "
            "$SLURM_TMPDIR/abm_opt if set, else ./output/workdirs. Node-local "
            "storage is much faster than the project filesystem."
        )
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

def small_scaffold_adjustment_aggrecan(value: float) -> float:
    """
    Adjusts the value for the small scaffold experimental data.
    The adjustment is based on the formula:
    value * (WORLD_VOLUME_MM3 / GEL_VOLUME_MM3)
    """
    return value * (WORLD_VOLUME_MM3 / GEL_VOLUME_MM3) # value in ug

EXPERIMENTAL_COLUMNS = {
    "cell_viability_percent": ("small_scaffold_cell_viability", None),
    "percent_diff": ("small_scaffold_percent_diff", None),
    "sGAG_total_ug": ("small_scaffold_aggrecan_ug", "aggrecan"),
    # Fold change in cell number relative to seeding; dimensionless, so no
    # conversion. Compared against the ABM's Total Cells at the metric's
    # day divided by Total Cells at its "relative_to_day". Present only at
    # 504 h (day 21) in experimental_config.csv.
    "fold_change": ("small_scaffold_fold_change", None),
}

def extract_small_scaffold_experimental(file_path: Path) -> pd.DataFrame:
    """
    Extract experimental data from the CSV, averaged per group and
    timepoint, with each measurement exposed under the small_scaffold_*
    name that OUTPUT_METRICS refers to.
    """
    df = pd.read_csv(file_path)
    df.columns = [c.strip() for c in df.columns]

    present = [c for c in EXPERIMENTAL_COLUMNS if c in df.columns]
    missing = [c for c in EXPERIMENTAL_COLUMNS if c not in df.columns]
    if missing:
        print(f"Note: experimental columns not in {file_path}, so any "
              f"OUTPUT_METRICS target using them will fail: {missing}")

    for col in present:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    # Average the replicates for each group and timepoint
    averages = df.groupby(["group", "time_hour"]).mean(numeric_only=True).reset_index()

    for col in present:
        target, conversion = EXPERIMENTAL_COLUMNS[col]
        if conversion == "aggrecan":
            averages[target] = averages[col].apply(small_scaffold_adjustment_aggrecan)
        elif conversion == "cells":
            averages[target] = averages[col].apply(small_scaffold_adjustment_cells)
        else:
            averages[target] = averages[col]

    return averages

def check_targets_against_experimental(
    experimental_df: pd.DataFrame, conditions: list[str]
) -> None:
    """
    Confirm every fitted target has a value at its day for every condition,
    before the optimizer starts. Without this a missing measurement shows
    up as a KeyError or a silent NaN many evaluations in.
    """
    problems = []
    for target in TARGETS:
        col = target["exp_column"]
        if col not in experimental_df.columns:
            problems.append(f"  {metric_label(target)}: no column '{col}' "
                            f"(see EXPERIMENTAL_COLUMNS)")
            continue
        for condition in conditions:
            key = (CONDITION_TO_GROUP[condition], target["day"] * 24)
            try:
                value = experimental_df.loc[key, col]
            except KeyError:
                problems.append(f"  {metric_label(target)}: no row for "
                                f"{CONDITION_TO_GROUP[condition]} at "
                                f"{target['day'] * 24} h (day {target['day']})")
                continue
            if pd.isna(value):
                problems.append(f"  {metric_label(target)}: '{col}' is blank "
                                f"for {CONDITION_TO_GROUP[condition]} at "
                                f"{target['day'] * 24} h")

    if problems:
        raise SystemExit(
            "Fitted targets have no experimental data to compare against:\n"
            + "\n".join(problems)
            + "\nEither add the measurements to experimental_config.csv, or "
              "drop 'exp_column' from those OUTPUT_METRICS entries to make "
              "them tracked-only."
        )
    print(f"All {len(TARGETS)} fitted target(s) have experimental values "
          f"for {len(conditions)} condition(s).")

def error(expected, actual):
    """
    Squared RELATIVE error between an experimental value and a simulated
    one -- the squared fractional difference, not the squared absolute
    difference:

        ((expected - actual) / expected) ** 2

    Relative error makes every target contribute on a comparable scale
    regardless of its units or magnitude. With absolute squared error the
    objective is dominated by whichever target has the largest numbers:
    percentages near 95 produce terms in the thousands, while a fold
    change near 1 or aggrecan near 0.07 ug produce terms near zero, so
    those targets end up effectively unfitted. In relative terms a 10%
    miss contributes 0.01 whichever target it came from.

    Consequences worth knowing:
      - total SSE now lands near 0.1-1 rather than in the thousands, so
        --tol must be scaled to match (default lowered to 1e-3).
      - a target whose expected value is legitimately near zero gets
        amplified; expected == 0 falls back to absolute error, since the
        fraction is undefined there.
    """
    if expected == 0:
        return (expected - actual) ** 2
    return ((expected - actual) / expected) ** 2

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

def _execute_one(job: dict) -> dict:
    """
    One ABM execution: write the config, run it, extract the metrics.

    Module-level and self-contained so it can be dispatched to a worker
    process. Each job gets its own working directory, so concurrent
    executions never share ./output/Output_Biomarkers.csv or a config file.
    """
    workdir = prepare_workdir(Path(job["workdir"]))
    config_path = workdir / "configFiles" / "simulation_config.json"

    create_sample_file(
        job["param_values"], Path(job["template"]), config_path,
        parameter_names=job["param_names"], condition=job["condition"],
    )
    biomarkers = run_ABM(config_path, numticks=NUMTICKS,
                         workdir=workdir, device=job["device"])

    return {
        "condition": job["condition"],
        "iteration": job["iteration"],
        "metrics": extract_output_metrics(biomarkers),
        "biomarkers": str(biomarkers),
    }

def run_conditions(
    conditions: list[str],
    experimental_df: pd.DataFrame,
    param_values: list,
    param_names: list,
    nfeval: int,
    args: argparse.Namespace,
) -> dict[str, dict[str, float]]:
    """
    Evaluate every condition for one parameter vector, running the
    (condition x --iters) ABM executions concurrently, and return the mean
    squared error per fitted target, per condition.

    """
    jobs = []
    for condition in conditions:
        for iteration in range(args.iters):
            jobs.append({
                "condition": condition,
                "iteration": iteration,
                "param_values": param_values,
                "param_names": param_names,
                "template": str(CONFIG_TEMPLATE.resolve()),
                "workdir": str(
                    args.workdir_root / f"eval{nfeval}_{condition}_{iteration}"
                ),
                "device": None,
                "numticks": args.numticks or NUMTICKS,
            })

    # Spread executions across the available GPUs round-robin, so they
    # don't all contend for device 0
    for i, job in enumerate(jobs):
        job["device"] = i % args.devices if args.devices > 1 else None

    workers = min(args.parallel, len(jobs))
    print(f"  Evaluation {nfeval}: {len(jobs)} ABM execution(s) "
          f"({len(conditions)} condition(s) x {args.iters} iter(s)), "
          f"{workers} at a time across {args.devices} GPU(s)")
    started = time.time()

    results = []
    if workers <= 1:
        for job in jobs:
            results.append(_execute_one(job))
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(_execute_one, job): job for job in jobs}
            for future in as_completed(futures):
                job = futures[future]
                try:
                    results.append(future.result())
                except Exception as exc:
                    # One failed execution invalidates this evaluation:
                    # averaging over the survivors would quietly change the
                    # objective's noise level mid-optimization
                    raise RuntimeError(
                        f"ABM execution failed (condition {job['condition']}, "
                        f"iteration {job['iteration']}, workdir "
                        f"{job['workdir']}): {exc}"
                    ) from exc

    # Collect errors per condition
    errors_by_condition: dict[str, dict[str, float]] = {}
    for condition in conditions:
        group_name = CONDITION_TO_GROUP[condition]
        runs = [r for r in results if r["condition"] == condition]

        per_target: dict[str, list[float]] = {metric_label(t): [] for t in TARGETS}
        for run in runs:
            simulated = run["metrics"]
            print(f"    {condition} iter {run['iteration']}: " + ", ".join(
                f"{label}={value:.4g}" for label, value in simulated.items()
            ))
            for target in TARGETS:
                label = metric_label(target)
                expected = float(
                    experimental_df.loc[
                        (group_name, target["day"] * 24), target["exp_column"]
                    ]
                )
                per_target[label].append(error(expected, simulated[label]))

        errors_by_condition[condition] = {
            label: float(np.mean(errs)) for label, errs in per_target.items()
        }

    elapsed = time.time() - started
    print(f"  Evaluation {nfeval}: {len(jobs)} execution(s) in "
          f"{elapsed / 60:.1f} min "
          f"({elapsed / len(jobs) / 60:.1f} min each, serialised)", flush=True)

    # Working directories are large (a 1009-row biomarker CSV plus logs per
    # execution) and there can be thousands of evaluations, so clear them
    # unless snapshots were asked for
    if not args.snapshots:
        for job in jobs:
            shutil.rmtree(job["workdir"], ignore_errors=True)

    return errors_by_condition

def spec_fingerprint() -> str:
    """
    Fingerprint of everything that changes what an SSE *means*: the fitted
    targets (their labels, days, columns and experimental columns) and the
    files the objective is computed from.

    This goes into the cache key, so editing OUTPUT_METRICS in ABM.py, the
    config template, or the experimental data invalidates old entries
    automatically instead of mixing incomparable SSEs -- or, worse,
    returning a cached dict that lacks the keys the new spec asks for.
    """
    payload = {
        "targets": [
            {
                "label": metric_label(t),
                "day": t["day"],
                "column": t["column"],
                "exp_column": t.get("exp_column"),
            }
            for t in TARGETS
        ],
        "numticks": NUMTICKS,
        "error_mode": "relative",
        "world_size_mm": list(WORLD_SIZE_MM),
        "gel_volume_mm3": GEL_VOLUME_MM3,
    }
    for f in (CONFIG_TEMPLATE, Path("experimental_config.csv")):
        try:
            payload[str(f)] = hashlib.sha256(Path(f).read_bytes()).hexdigest()[:16]
        except OSError:
            payload[str(f)] = None
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True).encode()
    ).hexdigest()[:16]


def cache_key(param_values: list, conditions: list[str], iters: int) -> str:
    """
    Cache identity for one objective evaluation. Includes the conditions
    and iteration count, since the same parameter vector evaluated over
    different conditions (or averaged over a different number of runs) is
    a different objective value -- and the spec fingerprint, so a changed
    OUTPUT_METRICS / template / experimental file can't produce a stale hit.

    Values are rounded to 10 significant figures so float formatting
    round-trips through JSON without missing a hit.
    """
    payload = json.dumps({
        "params": [float(f"{v:.10g}") for v in param_values],
        "conditions": conditions,
        "iters": iters,
        "spec": spec_fingerprint(),
    }, sort_keys=True)
    return hashlib.sha256(payload.encode()).hexdigest()[:32]


def load_cache(path: Path, enabled: bool) -> dict:
    if not enabled or not path.exists():
        return {}
    try:
        with open(path) as f:
            cache = json.load(f)
        print(f"Loaded {len(cache)} cached evaluation(s) from {path}")
        return cache
    except (json.JSONDecodeError, OSError) as exc:
        print(f"WARNING: could not read cache {path} ({exc}); starting fresh.")
        return {}


def save_cache(path: Path, cache: dict) -> None:
    """
    Write the cache atomically, so a job killed mid-write (e.g. by a
    walltime timeout) can't leave a truncated file that breaks the resume.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    with open(tmp, "w") as f:
        json.dump(cache, f)
    tmp.replace(path)

def make_objective(
    conditions: list[str],
    params: list[int],
    default_values: np.ndarray,
    param_names: list[str],
    experimental_df: pd.DataFrame,
    args: argparse.Namespace,
    cache: dict,
) -> Callable[[np.ndarray], float]:
    """
    Build the Nelder-Mead objective function.

    The returned function evaluates the ABM for every condition in
    `conditions` and returns the total SSE summed over all of them. Pass a
    single condition for a per-condition fit, or every condition for a
    joint fit -- the maths is the same, only the length of `conditions`
    differs.

    Evaluations are cached by parameter vector. Because Nelder-Mead is
    deterministic given the objective values, a resumed run replays the
    same sequence of points and gets cache hits until it reaches the point
    the previous job died on.
    """
state = {"nfeval": 1, "history": [], "cache_hits": 0, "started": time.time()}

    def objective(x: np.ndarray) -> float:
        nfeval = state["nfeval"]

        # Start from the spreadsheet defaults and overwrite only the
        # parameters being optimized with the simplex's current values
        sam = np.array(default_values, dtype=float)
        for idx, param in enumerate(params):
            sam[param] = x[idx]
        param_values = sam.tolist()

        key = cache_key(param_values, conditions, args.iters)
        cached = cache.get(key) if not args.no_cache else None

        # A hit from an older cache format (or one written before the
        # metric spec changed) may not carry every target. Treat that as a
        # miss rather than failing with a KeyError deep inside scipy.
        if cached is not None:
            needed = {metric_label(t) for t in TARGETS}
            have = cached.get("errors_by_condition", {})
            if any(needed - set(have.get(c, {})) for c in conditions):
                print(f"  Evaluation {nfeval}: stale cache entry "
                      f"(missing {sorted(needed - set(have.get(conditions[0], {})))}) "
                      f"-- re-running")
                cached = None

        if cached is not None:
            errors_by_condition = cached["errors_by_condition"]
            state["cache_hits"] += 1
            print(f"  Evaluation {nfeval}: cache hit "
                  f"(total SSE {cached['total_sse']:.6f}) -- skipping "
                  f"{len(conditions) * args.iters} ABM execution(s)")
        else:
            errors_by_condition = run_conditions(
                conditions, experimental_df, param_values,
                param_names, nfeval, args,
            )

        # Flatten to the objective vector Y: one entry per target per
        # condition, in (condition, target) order
        Y = np.array([
            errors_by_condition[condition][metric_label(target)]
            for condition in conditions
            for target in TARGETS
        ])

        total = float(np.sum(Y))

        if cached is None and not args.no_cache:
            cache[key] = {
                "errors_by_condition": errors_by_condition,
                "total_sse": total,
                "x": [float(v) for v in x],
                "spec": spec_fingerprint(),
            }
            # Written after every evaluation, not at the end: the whole
            # point is to survive the job being killed
            save_cache(args.cache, cache)

        # Progress line: one per evaluation, timestamped, with best-so-far
        # and a budget estimate. This is what you watch with
        # `tail -f slurm-<jobid>.out`.
        best = min([h["total_sse"] for h in state["history"]] + [total])
        elapsed = time.time() - state["started"]
        budget = args.maxfev or (len(params) + 1 + args.maxiter * 2)
        rate = elapsed / nfeval
        remaining = max(budget - nfeval, 0) * rate

        print(f"[{time.strftime('%H:%M:%S')}] eval {nfeval}/~{budget} "
              f"| SSE {total:.4g} | best {best:.4g} "
              f"| elapsed {elapsed / 60:.1f} min "
              f"| est. remaining {remaining / 60:.0f} min"
              + ("  (cached)" if cached is not None else ""), flush=True)

        print(formatted_string(nfeval, x, Y))
        print({c: errors_by_condition[c] for c in conditions})
        print(f"  total SSE over {len(conditions)} condition(s) "
              f"({', '.join(conditions)}): {total:.6f}", flush=True)

        state["history"].append({
            "nfeval": nfeval,
            "x": [float(v) for v in x],
            "errors_by_condition": errors_by_condition,
            "total_sse": total,
            "from_cache": cached is not None,
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
    cache: dict,
) -> dict:
    """
    Run one Nelder-Mead optimization fitting `conditions` together, and
    return a summary dict of the result.
    """
    label = "+".join(conditions)
    print(f"\n{'=' * 70}\nOptimizing conditions: {label}\n{'=' * 70}")

    objective = make_objective(
        conditions, params, default_values, param_names, experimental_df,
        args, cache,
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
    cache_hits = objective.state["cache_hits"]
    if cache_hits:
        print(f"  ({cache_hits}/{len(history)} evaluations replayed from cache)")
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
        "evaluations_from_cache": cache_hits,
        "success": bool(result.success),
        "message": str(result.message),
        "history": history,
    }


if __name__ == "__main__":
    args = parse_args()

    # Resolve the parallelism / storage defaults now, so everything
    # downstream sees concrete values
    if args.devices is None:
        args.devices = int(os.environ.get("SLURM_GPUS_ON_NODE", "1") or 1)
    if args.parallel is None:
        args.parallel = len(args.conditions) * args.iters
    if args.workdir_root is None:
        tmpdir = os.environ.get("SLURM_TMPDIR")
        args.workdir_root = (
            Path(tmpdir) / "abm_opt" if tmpdir else OUTPUT_DIR / "workdirs"
        )
    args.workdir_root.mkdir(parents=True, exist_ok=True)

    print(f"Parallelism: up to {args.parallel} concurrent ABM execution(s) "
          f"across {args.devices} GPU(s)")
    print(f"Working directories: {args.workdir_root}")
    print(f"World size: {WORLD_SIZE_MM[0]} x {WORLD_SIZE_MM[1]} x "
          f"{WORLD_SIZE_MM[2]} mm = {WORLD_VOLUME_MM3:.4g} mm^3 "
          f"({WORLD_VOLUME_MM3 / GEL_VOLUME_MM3:.3g} of the "
          f"{GEL_VOLUME_MM3:.0f} mm^3 gel)")
    print(f"Cache: {'disabled' if args.no_cache else args.cache}")

    cache = load_cache(args.cache, not args.no_cache)
    print(f"Objective spec fingerprint: {spec_fingerprint()} "
          f"(covers OUTPUT_METRICS, {CONFIG_TEMPLATE}, experimental_config.csv)")
    if cache:
        current = spec_fingerprint()
        stale = sum(1 for v in cache.values() if v.get("spec") not in (None, current))
        if stale:
            print(f"  {stale} cached entry/entries are from a different spec "
                  f"and will be ignored (safe to delete {args.cache}).")

    print(f"[{time.strftime('%H:%M:%S')}] Reading {PARAMETER_FILE} ...", flush=True)
    # Read parameter definitions
    df = pd.read_excel(PARAMETER_FILE)
    numpar = len(df)
    param_names = df["Parameter Name"].tolist()

    # Parameters to optimize: every row with a blank "Vary?" column
    params = extract_varying_param_indices(df)

    # Upfront cost estimate. Nelder-Mead's cost is driven by the number of
    # parameters, not by --maxiter: the initial simplex alone is one
    # evaluation per parameter plus one, before iteration 1 begins.
    n_par_opt = len(params)
    simplex_evals = n_par_opt + 1
    est_low = simplex_evals + args.maxiter
    est_high = simplex_evals + args.maxiter * 2
    n_exec = len(args.conditions) * args.iters
    if args.mode == "separate":
        n_exec = args.iters  # one condition at a time
    print(f"""
Cost estimate before starting:
  parameters optimized      : {n_par_opt}
  initial simplex           : {simplex_evals} evaluations (before iteration 1)
  --maxiter {args.maxiter:<4}              : ~{est_low}-{est_high} evaluations total"""
          + (f"\n  --maxfev {args.maxfev:<4}               : capped at {args.maxfev} evaluations"
             if args.maxfev else "")
          + f"""
  ABM executions per eval   : {n_exec} ({args.parallel} run concurrently)
  runs of the optimizer     : {1 if args.mode == 'joint' else len(args.conditions)} ({args.mode} mode)

  Multiply evaluations by the per-evaluation time printed below to
  estimate walltime. If that exceeds your allocation, the run will time
  out and resume from output/objective_cache.json on resubmission --
  but consider narrowing the "Vary?" column in parameters.xlsx, or
  --maxfev, instead.
""", flush=True)

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

    print(f"[{time.strftime('%H:%M:%S')}] Reading experimental_config.csv ...", flush=True)
    # Generate experimental values
    experimental_df = extract_small_scaffold_experimental(Path("experimental_config.csv"))
    experimental_df_indexed = experimental_df.set_index(['group', 'time_hour'])

    print("Experimental data:")
    print(experimental_df_indexed)

    check_targets_against_experimental(experimental_df_indexed, args.conditions)

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
            param_names, experimental_df_indexed, args, cache,
        )
        for group in condition_groups
    ]

    results = {
        "mode": args.mode,
        "parallel": args.parallel,
        "devices": args.devices,
        "conditions": list(args.conditions),
        "settings": {
            "iters_per_evaluation": args.iters,
            "maxiter": args.maxiter,
            "tol": args.tol,
            "numticks": NUMTICKS,
            "error_mode": "squared relative error",
            "parallel": args.parallel,
            "fitted_targets": FITTED_METRIC_KEYS,
            "tracked_metrics": [metric_label(m) for m in OUTPUT_METRICS],
        },
        "runs": runs,
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(results, f, indent=2)

    print(f"\n{'=' * 70}")
    print(f"[{time.strftime('%H:%M:%S')}] FINISHED. Mode: {args.mode}")
    for run in runs:
        print(f"  {'+'.join(run['conditions'])}: final SSE {run['final_sse']:.6f}"
              + (f" ({run['sse_percent_decrease']:.2f}% decrease)"
                 if run["sse_percent_decrease"] is not None else ""))
    if args.mode == "separate" and len(runs) > 1:
        print(f"  sum of per-condition SSEs: {sum(r['final_sse'] for r in runs):.6f}")
    print(f"Wrote {args.out}")
