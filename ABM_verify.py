import argparse
import json
from pathlib import Path
 
import numpy as np
import pandas as pd
 
from ABM import (
    create_sample_file,
    run_ABM,
    extract_output_metrics,
    get_template_parameter_paths,
    CONFIG_TEMPLATE,
    SAMPLE_CONFIG_PATH,
    OUTPUT_BIOMARKERS,
    OUTPUT_METRIC_KEYS,
    OUTPUT_METRICS,
    CONDITIONS,
    format_output_metrics,
)

# ======================
# Check that the ABM model is working as expected
# i.e. it produces expected outputs for known inputs
#
# Reads generated_samples_with_outputs.csv (written by
# ABM_generate_samples.py), re-runs the ABM for every row, and compares
# the fresh simulated outputs against the stored expected ones
#
# Sept 2026: Updated for JSON config structure
# ======================

# Columns written by ABM_generate_samples.py that are not parameters
METADATA_COLS = ["sample_id", "condition", "num_params_varied"]
 
# Legacy column names, kept so old CSVs are recognised as metadata rather
# than being mistaken for parameter columns
LEGACY_METADATA_COLS = ["param_set_id", "configuration", "num_params_mutated"]

def extract_param_names(
    df: pd.DataFrame,
    template_file: Path = CONFIG_TEMPLATE,
) -> list[str]:
    """
    Determine the parameter names, in the order create_sample_file expects
    them, and check that the CSV actually carries a column for each one.
 
    The template is the source of truth for ordering (same as
    parameters.xlsx, which create_sample_file cross-validates against)
    """
    template_paths = get_template_parameter_paths(template_file)
 
    missing = [p for p in template_paths if p not in df.columns]
    if missing:
        raise ValueError(
            f"generated_samples_with_outputs.csv is missing {len(missing)} "
            f"parameter column(s) present in '{template_file}':\n"
            + "\n".join(f"  {p}" for p in missing[:10])
            + "\nRegenerate the samples with the current ABM_generate_samples.py."
        )
 
    known = set(METADATA_COLS) | set(LEGACY_METADATA_COLS) | set(OUTPUT_METRIC_KEYS)
    unexpected = [
        c for c in df.columns
        if c not in known and c not in template_paths and not c.startswith("Unnamed")
    ]
    if unexpected:
        print(
            f"Note: ignoring {len(unexpected)} unrecognised column(s) in the CSV: "
            f"{unexpected[:10]}"
        )
 
    return template_paths

def resolve_condition(row: pd.Series) -> str | None:
    """
    Get the scaffold condition for a row. Returns None when the CSV
    predates the high/low split, in which case the template's own
    world_init values are left untouched
    """
    if "condition" not in row.index or pd.isna(row["condition"]):
        return None
 
    condition = str(row["condition"]).strip().lower()
    if condition not in CONDITIONS:
        raise ValueError(
            f"Unknown condition '{row['condition']}' in the samples CSV; "
            f"expected one of {CONDITIONS}."
        )
    return condition

def simulate_row(
    row: pd.Series,
    param_names: list[str],
    reps: int,
    template_file: Path,
) -> tuple[dict[str, float], list[dict[str, float]]]:
    """
    Rebuild this row's config, run the ABM `reps` times, and return the
    per-metric mean of the simulated outputs plus the individual runs
    """
    param_values = [float(row[name]) for name in param_names]
    condition = resolve_condition(row)
 
    runs = []
    for rep in range(reps):
        create_sample_file(
            param_values,
            template_file,
            SAMPLE_CONFIG_PATH,
            parameter_names=param_names,
            condition=condition,
        )
        run_ABM(SAMPLE_CONFIG_PATH)
        runs.append(extract_output_metrics(OUTPUT_BIOMARKERS))
 
    return {
        key: float(np.mean([run[key] for run in runs]))
        for key in OUTPUT_METRIC_KEYS
    }, runs

def summarise(errors_by_metric: dict[str, list[float]],
              expected_by_metric: dict[str, list[float]]) -> dict:
    """
    Build the summary metrics block: MSE, RMSE, MAE and normalised RMSE
    (RMSE as a fraction of the mean expected value) per output variable.
    NRMSE is the useful one for comparing across metrics, since collagen
    (ug) and cell counts differ by orders of magnitude
    """
    summary = {"mean_squared_error": {}, "rmse": {}, "mean_absolute_error": {}, "nrmse": {}}
 
    for key in OUTPUT_METRIC_KEYS:
        errs = np.asarray(errors_by_metric[key], dtype=float)
        expected = np.asarray(expected_by_metric[key], dtype=float)
 
        mse = float(np.mean(errs ** 2))
        rmse = float(np.sqrt(mse))
        expected_mean = float(np.mean(np.abs(expected)))
 
        summary["mean_squared_error"][key] = mse
        summary["rmse"][key] = rmse
        summary["mean_absolute_error"][key] = float(np.mean(np.abs(errs)))
        summary["nrmse"][key] = rmse / expected_mean if expected_mean else None
 
    return summary

def verify(
    samples_csv: Path,
    output_json: Path,
    reps: int,
    template_file: Path,
    limit: int | None = None,
) -> dict:
    samples_df = pd.read_csv(samples_csv)
    if limit is not None:
        samples_df = samples_df.head(limit)
 
    param_names = extract_param_names(samples_df, template_file)
 
    missing_metrics = [k for k in OUTPUT_METRIC_KEYS if k not in samples_df.columns]
    if missing_metrics:
        raise ValueError(
            f"Samples CSV has no expected-output column(s) for {missing_metrics}. "
            f"extract_output_metrics in ABM.py and the CSV are out of sync -- "
            f"regenerate the samples."
        )
 
    print(
        f"Verifying {len(samples_df)} sample(s) from {samples_csv} "
        f"with {len(param_names)} parameters, {reps} rep(s) each."
    )
 
    per_sample_results = []
    errors_by_metric = {key: [] for key in OUTPUT_METRIC_KEYS}
    expected_by_metric = {key: [] for key in OUTPUT_METRIC_KEYS}
 
    for idx, row in samples_df.iterrows():
        simulated, runs = simulate_row(row, param_names, reps, template_file)
 
        errors, squared_errors, relative_errors = {}, {}, {}
        for key in OUTPUT_METRIC_KEYS:
            expected = float(row[key])
            sim = simulated[key]
            err = sim - expected
 
            errors[key] = err
            squared_errors[key] = err ** 2
            relative_errors[key] = err / expected if expected else None
 
            errors_by_metric[key].append(err)
            expected_by_metric[key].append(expected)
 
        result = {
            "row_index": int(idx),
            "sample_id": int(row["sample_id"]) if "sample_id" in row.index else int(idx),
            "condition": resolve_condition(row),
            "input_parameters": {name: float(row[name]) for name in param_names},
            "expected_outputs": {k: float(row[k]) for k in OUTPUT_METRIC_KEYS},
            "simulated_outputs": simulated,
            "errors": errors,
            "squared_errors": squared_errors,
            "relative_errors": relative_errors,
        }
        if reps > 1:
            result["replicate_outputs"] = runs
 
        per_sample_results.append(result)
 
        print(
            f"  row {idx} (sample {result['sample_id']}, {result['condition']}): "
            + ", ".join(
                f"{k}: sim={simulated[k]:.4g} vs exp={float(row[k]):.4g}"
                for k in OUTPUT_METRIC_KEYS
            )
        )
 
    summary = summarise(errors_by_metric, expected_by_metric)
    summary["total_samples"] = len(per_sample_results)
    summary["reps_per_sample"] = reps
 
    results_json = {
        "config": {
            "samples_csv": str(samples_csv),
            "config_template": str(template_file),
            "output_metrics": [
                {"label": k, "day": m["day"], "column": m["column"]}
                for k, m in zip(OUTPUT_METRIC_KEYS, OUTPUT_METRICS)
            ],
            "num_parameters": len(param_names),
        },
        "per_sample_results": per_sample_results,
        "summary_metrics": summary,
    }
 
    with open(output_json, "w") as f:
        json.dump(results_json, f, indent=2)
 
    print(f"\nWrote {output_json}")
    print("NRMSE by output (lower is better):")
    for key in OUTPUT_METRIC_KEYS:
        nrmse = summary["nrmse"][key]
        print(f"  {key}: {'n/a' if nrmse is None else f'{nrmse:.4f}'}")
 
    return results_json

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Verify the ABM reproduces known outputs for known inputs"
    )
    parser.add_argument(
        "--samples", type=Path, default=Path("generated_samples_with_outputs.csv"),
        help="CSV of parameter sets and expected outputs from ABM_generate_samples.py"
    )
    parser.add_argument(
        "--out", type=Path, default=Path("verification_results.json"),
        help="Where to write the verification results"
    )
    parser.add_argument(
        "--reps", type=int, default=1,
        help="ABM runs per sample, averaged before scoring (the ABM is stochastic)"
    )
    parser.add_argument(
        "--template", type=Path, default=CONFIG_TEMPLATE,
        help="JSON config template to build each sample's config from"
    )
    parser.add_argument(
        "--limit", type=int, default=None,
        help="Only verify the first N rows (useful for a quick test)"
    )
    args = parser.parse_args()
 
    verify(
        samples_csv=args.samples,
        output_json=args.out,
        reps=args.reps,
        template_file=args.template,
        limit=args.limit,
    )