import argparse
import json
from pathlib import Path

from ABM import (
    create_sample_file,
    run_ABM,
    extract_output_metrics,
    CONFIG_TEMPLATE,
    SAMPLE_CONFIG_PATH,
    OUTPUT_BIOMARKERS,
    CONDITIONS,
    PARAMETER_FILE,
    format_output_metrics,
    metric_label,
    OUTPUT_METRICS,
)

import numpy as np
import pandas as pd
from scipy.stats import truncnorm

# ==================================================
# Script to generate samples, including parameter sets
# and expected outputs, for the ABM model.
#
# Sept 2026: Updated to work with the JSON config structure.
# A single simulation_config.template.json (with mutable params
# under "biology") is used to write a merged JSON config file per sample
# instead of the old Sample.txt
#
# Split into phases so the ABM runs can be split up across cluster
# jobs (see ABM_generate_samples_job.sh)
#
#   write-configs  sample the parameters and write one JSON config per
#                  (sample, condition) task, plus sample_parameters.csv
#   run --task N   run the ABM for task N only and write its metrics to
#                  output/task_metrics/; one cluster job per task
#   collect        join sample_parameters.csv with every task's metrics
#                  into generated_samples_with_outputs.csv
#   all            do everything serially in one process (the old
#                  behaviour); fine for a couple of samples locally
#
# ==================================================

# --- Settings to edit -------------------------------------------------
input_file = PARAMETER_FILE
sheet_name = "Sheet1"
n_par = 67  # Number of parameters
rng_seed = None  # Set an int here for reproducibility, or leave None

SAMPLES_CONFIG_DIR = Path("configFiles/samples")       # written by write-configs
SAMPLE_PARAMETERS_CSV = Path("sample_parameters.csv")  # written by write-configs
TASK_METRICS_DIR = Path("output/task_metrics")         # written by run
SAMPLES_OUTPUT_CSV = Path("generated_samples_with_outputs.csv")  # written by collect
# -----------------------------------------------------------------------

rng = np.random.default_rng(rng_seed)

def sample_truncnorm(mean, sigma, lower, upper, size, rng):
    """
    Draw `size` samples from a normal(mean, sigma) distribution truncated
    to [lower, upper], so that bounds [a, b] are a hard min/max.
    """
    a_std = (lower - mean) / sigma
    b_std = (upper - mean) / sigma
    return truncnorm.rvs(a_std, b_std, loc=mean, scale=sigma, size=size, random_state=rng)

def generate_param_df() -> pd.DataFrame:
    """
    Generate a DataFrame of parameters with their bounds and default values from spreadsheet.
    
    Columns:
        'Parameter Number': 1-indexed row/param number
        'Parameter Name': colon-joined path matching the config template JSON
        'Mean': default/mean value; used directly when Vary? is 'N'
        'Lower': lower bound (a) of that param
        'Upper': upper bound (b) of that param
        'SD': std dev, precomputed as (Upper - Lower)/2
        'Type': distribution type: 1=lognormal, 2=discrete lognormal, 3=normal, 4=uniform 
        'Vary?': blank=vary this parameter, 'N'=use Mean for every run
    """
    df = pd.read_excel(input_file, sheet_name=sheet_name)
    df = df.iloc[:n_par].reset_index(drop=True)
    return df

def mutate_parameters(params: pd.DataFrame) -> pd.DataFrame:
    """
    Mutate the given parameters in the DataFrame.
    Rows marked 'N' in the 'Vary?' column are held at their Mean value
    instead of randomly sampled
    """
    mutated_params = params.copy()
    for i in range(len(mutated_params)):
        row = mutated_params.iloc[i]
        a = row["Lower"]
        b = row["Upper"]
        sd = row["SD"]
        dist_type = row["Type"]
        mean_val = row["Mean"]
        
        vary_raw = row["Vary?"]
        dont_vary = isinstance(vary_raw, str) and vary_raw.strip().upper() == "N"
        
        if dont_vary:
            mutated_params.loc[i, "value"] = mean_val
            continue
        
        if dist_type == 1:
            # lognormal, truncated to [a, b]
            mu = np.log(mean_val)
            sigma = (np.log(b) - np.log(a)) / 2
            log_r = sample_truncnorm(mu, sigma, np.log(a), np.log(b), 1, rng)
            mutated_params.loc[i, "value"] = np.exp(log_r)[0]
        elif dist_type == 2:
            # discrete lognormal, rounded to nearest 0.5
            mu = np.log(mean_val)
            sigma = (np.log(b) - np.log(a)) / 2
            log_x = sample_truncnorm(mu, sigma, np.log(a), np.log(b), 1, rng)
            x = np.exp(log_x)[0]
            y = np.floor(x)
            frac = x - y
            r = y + 0.5 if frac >= 0.5 else y
            mutated_params.loc[i, "value"] = np.clip(r, a, b)
        elif dist_type == 3:
            # normal, truncated to [a, b]
            r = sample_truncnorm(mean_val, sd, a, b, 1, rng)
            mutated_params.loc[i, "value"] = r[0]
        elif dist_type == 4:
            # uniform between lower and upper bound
            mutated_params.loc[i, "value"] = rng.uniform(a, b)
        else:
            raise ValueError(f"Unknown distribution type '{dist_type}' at parameter {i + 1}")
            
    return mutated_params

def task_list(n_samples: int) -> list[tuple[int, str]]:
    """
    The canonical (sample_id, condition) ordering. Task N is
    task_list(n)[N], so the cluster array index maps to exactly one ABM
    run. Every phase derives its task numbering from this one function --
    don't reorder it between phases of the same sweep.
    """
    return [(i, condition) for i in range(n_samples) for condition in CONDITIONS]

def task_name(sample_id: int, condition: str) -> str:
    """Stable filename stem for a task, e.g. 'sample0_high'."""
    return f"sample{sample_id}_{condition}"

def write_configs(n_samples: int) -> pd.DataFrame:
    """
    Phase 1: draw parameter sets and write one ready-to-run JSON config
    per task, plus sample_parameters.csv recording the drawn values
 
    Each sample's parameters are drawn ONCE and shared by all of its
    conditions, so the high/low runs of a given sample differ only in the
    scaffold fields
    """
    param_df = generate_param_df()
    param_names = param_df["Parameter Name"].tolist()
    num_varying = param_df["Vary?"].apply(
        lambda v: not (isinstance(v, str) and v.strip().upper() == "N")).sum()
 
    SAMPLES_CONFIG_DIR.mkdir(parents=True, exist_ok=True)
 
    rows = []
    for i in range(n_samples):
        sample_df = mutate_parameters(param_df)
        param_values = sample_df["value"].tolist()
 
        for condition in CONDITIONS:
            config_path = SAMPLES_CONFIG_DIR / f"{task_name(i, condition)}.json"
            create_sample_file(
                param_values,
                CONFIG_TEMPLATE,
                config_path,
                parameter_names=param_names,
                condition=condition,
            )
            rows.append({
                "task": len(rows),
                "sample_id": i,
                "condition": condition,
                "config_file": str(config_path),
                "num_params_varied": num_varying,
                **{name: param_values[j] for j, name in enumerate(param_names)},
            })
 
    df = pd.DataFrame(rows)
    df.to_csv(SAMPLE_PARAMETERS_CSV, index=False)
 
    expected = task_list(n_samples)
    assert [(r["sample_id"], r["condition"]) for r in rows] == expected, \
        "config ordering does not match task_list()"
 
    print(f"Wrote {len(rows)} config(s) to {SAMPLES_CONFIG_DIR} "
          f"and {SAMPLE_PARAMETERS_CSV}")
    print(f"Task indices 0..{len(rows) - 1}")
    return df

def run_task(task: int) -> dict[str, float]:
    """
    Phase 2: run the ABM for a single task and write its metrics to
    output/task_metrics/<task_name>.json
 
    Reads the config written by write-configs. Intended to be
    the body of one cluster job; run from a private working directory so
    concurrent tasks don't overwrite each other's Output_Biomarkers.csv.
    """
    params_df = pd.read_csv(SAMPLE_PARAMETERS_CSV)
    matching = params_df[params_df["task"] == task]
    if matching.empty:
        raise ValueError(
            f"No task {task} in {SAMPLE_PARAMETERS_CSV} "
            f"(valid range 0..{len(params_df) - 1}). "
            f"Run 'write-configs' first, or check the array indices."
        )
    row = matching.iloc[0]
 
    config_path = Path(row["config_file"])
    if not config_path.exists():
        raise FileNotFoundError(
            f"Config {config_path} for task {task} is missing; re-run "
            f"'write-configs' (the configs and CSV must come from the same sweep)."
        )
 
    print(f"Task {task}: sample {row['sample_id']}, condition "
          f"{row['condition']}, config {config_path}")
 
    run_ABM(config_path)
    metrics = extract_output_metrics(OUTPUT_BIOMARKERS)
    print(f"Task {task}: {format_output_metrics(metrics)}")
 
    TASK_METRICS_DIR.mkdir(parents=True, exist_ok=True)
    metrics_path = TASK_METRICS_DIR / f"{task_name(int(row['sample_id']), row['condition'])}.json"
    with open(metrics_path, "w") as f:
        json.dump({
            "task": task,
            "sample_id": int(row["sample_id"]),
            "condition": row["condition"],
            "metrics": metrics,
        }, f, indent=2)
 
    print(f"Task {task}: wrote {metrics_path}")
    return metrics

def collect(strict: bool = False) -> pd.DataFrame:
    """
    Phase 3: join sample_parameters.csv with each task's metrics JSON into
    generated_samples_with_outputs.csv, in the same column layout the
    serial version produced
 
    Tasks with no metrics file are reported and skipped. 
    Pass strict=True to fail instead if you intend to verify against 
    the full set
    """
    params_df = pd.read_csv(SAMPLE_PARAMETERS_CSV)
    metric_keys = [metric_label(m) for m in OUTPUT_METRICS]
 
    rows, missing = [], []
    for _, row in params_df.iterrows():
        stem = task_name(int(row["sample_id"]), row["condition"])
        metrics_path = TASK_METRICS_DIR / f"{stem}.json"
 
        if not metrics_path.exists():
            missing.append((int(row["task"]), stem))
            continue
 
        with open(metrics_path) as f:
            metrics = json.load(f)["metrics"]
 
        absent = [k for k in metric_keys if k not in metrics]
        if absent:
            raise ValueError(
                f"{metrics_path} is missing metric(s) {absent}. It was "
                f"probably written before OUTPUT_METRICS changed -- re-run "
                f"the affected tasks."
            )
 
        # Drop the phase-plumbing columns so the output matches what the
        # serial path produced
        record = row.drop(labels=["task", "config_file"]).to_dict()
        rows.append({**record, **{k: metrics[k] for k in metric_keys}})
 
    if missing:
        print(f"WARNING: {len(missing)} task(s) have no metrics file:")
        for task, stem in missing[:20]:
            print(f"  task {task} ({stem})")
        if strict:
            raise SystemExit(
                f"{len(missing)} task(s) incomplete and --strict was given."
            )
 
    if not rows:
        raise SystemExit("No completed tasks found -- nothing to collect.")
 
    result_df = pd.DataFrame(rows)
    result_df.to_csv(SAMPLES_OUTPUT_CSV, index=False)
    print(f"Collected {len(rows)}/{len(params_df)} task(s) into {SAMPLES_OUTPUT_CSV}")
    return result_df

def generate_samples(n_samples: int) -> pd.DataFrame:
    """
    Run all three phases serially in one process (the original
    behaviour) (fine for a few samples locally or something)
    """
    write_configs(n_samples)
    for task, _ in enumerate(task_list(n_samples)):
        run_task(task)
    return collect()

def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Generate ABM parameter samples and their expected outputs."
    )
    sub = parser.add_subparsers(dest="phase", required=True)
 
    p_write = sub.add_parser(
        "write-configs",
        help="Draw parameter sets and write one JSON config per task (no ABM runs)."
    )
    p_write.add_argument("--n-samples", type=int, default=10,
                         help="Number of parameter sets (each run once per condition).")
 
    p_run = sub.add_parser("run", help="Run the ABM for a single task.")
    p_run.add_argument("--task", type=int, required=True,
                       help="0-based task index (see write-configs output).")
 
    p_collect = sub.add_parser(
        "collect", help="Join the per-task metrics into the samples CSV."
    )
    p_collect.add_argument("--strict", action="store_true",
                           help="Fail if any task is missing its metrics file.")
 
    p_all = sub.add_parser("all", help="Run every phase serially in this process.")
    p_all.add_argument("--n-samples", type=int, default=10)
 
    p_count = sub.add_parser(
        "count", help="Print the number of tasks for a given sample count and exit."
    )
    p_count.add_argument("--n-samples", type=int, default=10)
 
    args = parser.parse_args(argv)
 
    if args.phase == "write-configs":
        write_configs(args.n_samples)
    elif args.phase == "run":
        run_task(args.task)
    elif args.phase == "collect":
        collect(strict=args.strict)
    elif args.phase == "all":
        generate_samples(args.n_samples)
    elif args.phase == "count":
        print(len(task_list(args.n_samples)))
 
 
if __name__ == "__main__":
    main()