import os
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

# Regexes for logated the variable numeric fields
# under the "biology" key in the JSON config template
# for ABMs
KEY_OPEN_RE = re.compile(r'^\s*"(?P<key>[^"]+)"\s*:\s*\{\s*$')
CLOSE_RE = re.compile(r'^\s*\}')
LINE_RE = re.compile(
    r'^(?P<prefix>\s*"(?P<key>[^"]+)"\s*:\s*)'
    r'(?P<value>-?\d+\.?\d*)'
    r'(?P<suffix>,?\s*//.*)$'
)

# Looks for direct fields (e.g. world_init scaffold params in the JSON
# that aren't tagged with // param identifiers like the biology params)
LEAF_ANY_RE = re.compile(
    r'^(?P<prefix>\s*"(?P<key>[^"]+)"\s*:\s*)'
    r'(?P<value>-?\d+\.?\d*)'
    r'(?P<suffix>,?\s*(?://.*)?)$'
)

# Scaffold condition: sets world_init.alginate.high_mw_ratio / low_mw_ratio
# to switch between the 'high' and 'low' MW scaffold conditions
CONDITION_FIELDS = {
    "high": {"world_init:alginate:high_mw_ratio": 1, "world_init:alginate:low_mw_ratio": 0},
    "low": {"world_init:alginate:high_mw_ratio": 0, "world_init:alginate:low_mw_ratio": 1},
}

# ---------------------------------------------------------------------
# Shared paths/constants
#
# Every script in this repo (ABM_generate_samples, ABM_optimize,
# ABM_verify) imports these when necessary
# ---------------------------------------------------------------------

CONFIG_TEMPLATE = Path("simulation_config.template.json")
SAMPLE_CONFIG_PATH = Path("configFiles/simulation_config_sample.json")
OUTPUT_DIR = Path("output")
OUTPUT_BIOMARKERS = OUTPUT_DIR / "Output_Biomarkers.csv"
PARAMETER_FILE = Path("parameters.xlsx")
CONDITIONS = list(CONDITION_FIELDS.keys())

CLOCK_COL = "clock (30 min)"
TICKS_PER_DAY = 48

# Column names taken directly from the Output_Biomarkers.csv header row
CLOCK_COL = "clock (30 min)"
COLLAGEN_COL = "Collagen (ug)"
TOTAL_CELLS_COL = "Total Cells"
AGGRECAN_COL = "Aggrecan (ug)"
CELL_VIABILITY_COL = "Viability Rate(%)"
PERCENT_DIFF_COL = "Differentiation (%)"

# 30-minute ticks, so 48 ticks per simulated day
TICKS_PER_DAY = 48

# ---------------------------------------------------------------------
# This list of output metrics is shared across all the files in the
# optimization protocol. ONLY edit here whenever you want to change what is
# read out of the ABM!

OUTPUT_METRICS = [
    # Day 7 -- fitted against experimental data
    {"day": 7, "column": CELL_VIABILITY_COL, "label": "day_7_cell_viability",
     "exp_column": "small_scaffold_cell_viability"},
    {"day": 7, "column": PERCENT_DIFF_COL, "label": "day_7_percent_diff",
     "exp_column": "small_scaffold_percent_diff"},
    # Day 7 -- tracked only (no experimental data at this timepoint)
    {"day": 7, "column": AGGRECAN_COL, "label": "day_7_aggrecan"},
    {"day": 7, "column": TOTAL_CELLS_COL, "label": "day_7_cells"},

    # Day 21 -- fitted against experimental data
    {"day": 21, "column": AGGRECAN_COL, "label": "day_21_aggrecan",
     "exp_column": "small_scaffold_aggrecan_ug"},
    {"day": 21, "column": CELL_VIABILITY_COL, "label": "day_21_cell_viability",
     "exp_column": "small_scaffold_cell_viability"},
    {"day": 21, "column": PERCENT_DIFF_COL, "label": "day_21_percent_diff",
     "exp_column": "small_scaffold_percent_diff"},
    # Day 21 -- tracked only
    {"day": 21, "column": TOTAL_CELLS_COL, "label": "day_21_cells"},
]

# ---------------------------------------------------------------------

def _slug(column: str) -> str:
    """
    Turn a biomarker column name into a label-safe slug, e.g.
    "Collagen (ug)" -> "collagen", "Viability Rate(%)" -> "viability_rate".
    """
    base = re.sub(r"\(.*?\)", "", column)  # drop units
    base = re.sub(r"[^0-9A-Za-z]+", "_", base).strip("_")
    return base.lower()

def metric_label(metric: dict) -> str:
    """
    The label for a metric spec: its explicit "label" if given, otherwise
    "day_<day>_<column slug>".
    """
    return metric.get("label") or f"day_{metric['day']}_{_slug(metric['column'])}"

# Metric labels, in order and derived from OUTPUT_METRICS
OUTPUT_METRIC_KEYS = [metric_label(m) for m in OUTPUT_METRICS]

# subset of metrics that ABM_optimize fits: those with experimental data to
# compare against. Each contributes one squared-error term per condition
FITTED_METRICS = [m for m in OUTPUT_METRICS if m.get("exp_column")]
FITTED_METRIC_KEYS = [metric_label(m) for m in FITTED_METRICS]

# ticks needed to reach the latest day in OUTPUT_METRICS (+1 because tick
# numbering starts at 0)
METRICS_NUMTICKS = TICKS_PER_DAY * max(m["day"] for m in OUTPUT_METRICS) + 1


def load_template_lines(path: Path) -> list[str]:
    with open(path, "r") as f:
        return f.readlines()

def find_variable_line_entries(lines: list[str]) -> list[tuple[int, str]]:
    """
    Walk the template's nested JSON object keys (via brace tracking) and
    collect every //-tagged, mutable numeric field found under the
    "biology" key

    Returns a list of (line_index, parameter_path) tuples, where
    parameter_path is a colon-joined path of the field's ancestor keys
    below "biology" plus its own key -- e.g. for:

        "biology": {
          "cell": {
            "proliferation": {
              "hours_between_proliferation": 24, // k2

    the entry would be (line_index, "cell:proliferation:hours_between_proliferation").

    This path is meant to match the "Parameter Name" column in
    parameters.xlsx, so the two can be cross-validated
    """
    entries = []
    stack: list[str] = []
    in_biology = False

    for i, line in enumerate(lines):
        m_open = KEY_OPEN_RE.match(line)
        if m_open:
            key = m_open.group("key")
            stack.append(key)
            if key == "biology":
                in_biology = True
            continue

        if CLOSE_RE.match(line):
            if stack:
                popped = stack.pop()
                if popped == "biology":
                    in_biology = False
            continue

        if in_biology:
            m_val = LINE_RE.match(line)
            if m_val:
                path = ":".join(stack[1:] + [m_val.group("key")])  # drop leading "biology"
                entries.append((i, path))

    return entries

def find_variable_line_indices(lines: list[str]) -> list[int]:
    """
    Backwards-compatible helper returning just the line indices (no paths)
    """
    return [idx for idx, _ in find_variable_line_entries(lines)]

def get_template_parameter_paths(template_file: Path) -> list[str]:
    """
    Return the ordered list of parameter_path values for a template's
    //-tagged "biology" fields. Useful for generating or validating the
    "Parameter Name" column of parameters.xlsx
    """
    lines = load_template_lines(template_file)
    return [path for _, path in find_variable_line_entries(lines)]

def find_all_leaf_entries(lines: list[str]) -> list[tuple[int, str]]:
    """
    Like find_variable_line_entries, but walks the ENTIRE template (not
    just under "biology") and matches any numeric leaf field, whether or
    not it has a trailing // comment. Returns (line_index, full_path)
    tuples, where full_path includes every ancestor key from the root of
    the file -- e.g. "world_init:alginate:high_mw_ratio".
    """
    entries = []
    stack: list[str] = []

    for i, line in enumerate(lines):
        m_open = KEY_OPEN_RE.match(line)
        if m_open:
            stack.append(m_open.group("key"))
            continue

        if CLOSE_RE.match(line):
            if stack:
                stack.pop()
            continue

        m_val = LEAF_ANY_RE.match(line)
        if m_val:
            path = ":".join(stack + [m_val.group("key")])
            entries.append((i, path))

    return entries

def set_leaf_value(lines: list[str], full_path: str, value: float) -> list[str]:
    """
    Set a single numeric leaf field, identified by its full colon-joined
    path from the root of the file (e.g.
    "world_init:alginate:high_mw_ratio"), to the given value. Returns a
    new list of lines with that one field updated; raises ValueError if
    the path isn't found or is ambiguous.
    """
    entries = find_all_leaf_entries(lines)
    matches = [idx for idx, path in entries if path == full_path]

    if not matches:
        raise ValueError(f"Could not find field '{full_path}' in template.")
    if len(matches) > 1:
        raise ValueError(f"Field '{full_path}' matched multiple lines: {matches}")

    idx = matches[0]
    new_lines = lines.copy()
    m = LEAF_ANY_RE.match(new_lines[idx])
    prefix, suffix = m.group("prefix"), m.group("suffix") or ""
    val_str = str(int(value)) if float(value).is_integer() else f"{value:.6g}"
    new_lines[idx] = f"{prefix}{val_str}{suffix}\n"
    return new_lines

def apply_scaffold_condition(lines: list[str], condition: str) -> list[str]:
    """
    Apply the "high" or "low" molecular-weight scaffold condition by
    setting world_init.alginate.high_mw_ratio / low_mw_ratio to (1, 0)
    or (0, 1) respectively.
    """
    if condition not in CONDITION_FIELDS:
        raise ValueError(f"Unknown condition '{condition}'; expected one of {list(CONDITION_FIELDS)}")

    new_lines = lines
    for path, value in CONDITION_FIELDS[condition].items():
        new_lines = set_leaf_value(new_lines, path, value)
    return new_lines

def create_sample_file(
    parameter_list: list[float],
    template_file: Path,
    out_path: Path,
    parameter_names: list[str] | None = None,
    condition: str | None = None,
    ) -> None:
    """
    Create a sample file (JSON template) based on the given parameters and
    the condition (high/low)
    All other parts of the template (the rest of world_init, chemistry) are
    copied through unchanged
    """
    lines = load_template_lines(template_file)
    entries = find_variable_line_entries(lines)
    var_idxs = [idx for idx, _ in entries]
    var_paths = [path for _, path in entries]

    if len(parameter_list) != len(var_idxs):
        raise ValueError(
            f"Mismatch: got {len(parameter_list)} parameter values but "
            f"template '{template_file}' has {len(var_idxs)} //-tagged "
            f"parameters under 'biology'. Check that the parameter order/count "
            f"matches the template's tagged-line order."
        )

    if parameter_names is not None:
        if len(parameter_names) != len(var_paths):
            raise ValueError(
                f"Mismatch: got {len(parameter_names)} parameter names but "
                f"template '{template_file}' has {len(var_paths)} //-tagged parameters."
            )
        mismatches = [
            (i, given, expected)
            for i, (given, expected) in enumerate(zip(parameter_names, var_paths))
            if given != expected
        ]
        if mismatches:
            details = "\n".join(
                f"  row {i}: parameters.xlsx has '{given}', template has '{expected}'"
                for i, given, expected in mismatches[:10]
            )
            raise ValueError(
                f"Parameter Name mismatch between parameters.xlsx and template "
                f"'{template_file}' ({len(mismatches)} mismatched row(s)):\n{details}\n"
                f"Check that parameters.xlsx rows are in the same order as the "
                f"template's //-tagged lines."
            )

    new_lines = lines.copy()
    for idx, val in zip(var_idxs, parameter_list):
        m = LINE_RE.match(new_lines[idx])
        prefix, suffix = m.group("prefix"), m.group("suffix") or ""
        val_str = str(int(val)) if float(val).is_integer() else f"{val:.6g}"
        new_lines[idx] = f"{prefix}{val_str}{suffix}\n"

    if condition is not None:
        new_lines = apply_scaffold_condition(new_lines, condition)

    if out_path.exists():
        os.remove(out_path)

    with open(out_path, "w") as f:
        f.writelines(new_lines)

def extract_biomarkers_at_day(df: pd.DataFrame, day: int, source: str = "") -> pd.Series:
    """
    Return the Output_Biomarkers row for a given simulated day, looked up
    by the clock column
    """
    tick = TICKS_PER_DAY * day
    match = df[df[CLOCK_COL] == tick]
    if match.empty:
        raise ValueError(
            f"No row in {source or 'the biomarker file'} with {CLOCK_COL} == "
            f"{tick} (day {day}). Check that run_ABM's --numticks covers day "
            f"{day} and that the simulation completed -- see output/stderr.txt."
        )
    return match.iloc[0]


def extract_output_metrics(output_file: Path,
                           metrics: list[dict] | None = None,
                           ) -> dict[str, float]:
    """
    Extract output metrics from the ABM output_biomarkers.csv file.

    Reads the (day, column) pairs declared in OUTPUT_METRICS at the top of
    this file. Don't edit this function! If you need to change what
    metrics the pipeline scores, edit OUTPUT_METRICS directly!

    Pass 'metrics' to override for a one-off extraction of a new metric
    that you didn't declare in OUTPUT_METRICS.
    """
    metrics = OUTPUT_METRICS if metrics is None else metrics

    df = pd.read_csv(output_file)
    df.columns = [c.strip() for c in df.columns]

    missing = sorted({m["column"] for m in metrics} - set(df.columns))
    if missing:
        raise ValueError(
            f"OUTPUT_METRICS refers to column(s) not present in {output_file}: "
            f"{missing}.\nAvailable columns: {list(df.columns)}"
        )

    # Fetch each needed day once, then pull every column for that day
    rows = {
        day: extract_biomarkers_at_day(df, day, str(output_file))
        for day in sorted({m["day"] for m in metrics})
    }

    return {
        metric_label(m): float(rows[m["day"]][m["column"]])
        for m in metrics
    }

def format_output_metrics(values: dict[str, float]) -> str:
    """
    One-line, human-readable rendering of an extract_output_metrics result
    """
    return ", ".join(f"{label}={values[label]:.4g}" for label in values)

def prepare_workdir(workdir: Path, extra_files: list[Path] | None = None) -> Path:
    """
    Set up a private working directory for one ABM execution (helps with concurrent
    executions).
    """
    workdir = Path(workdir)
    (workdir / "output").mkdir(parents=True, exist_ok=True)
    (workdir / "configFiles").mkdir(parents=True, exist_ok=True)

    bin_link = workdir / "bin"
    if not bin_link.exists():
        bin_link.symlink_to(Path("bin").resolve(), target_is_directory=True)

    for f in extra_files or []:
        target = workdir / Path(f).name
        if not target.exists():
            shutil.copy(Path(f), target)

    return workdir

def run_ABM(config_file: Path, numticks: int = METRICS_NUMTICKS, workdir: Path | None = None, device: int | None = None,) -> Path:
    """
    Run the ABM with a given JSON configuration file.
    """
    base = Path(workdir) if workdir is not None else Path(".")
    out_dir = base / OUTPUT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)

    stdout_file_name = OUTPUT_DIR / "stdout.txt"
    stderr_file_name = OUTPUT_DIR / "stderr.txt"
    biomarkers = base / OUTPUT_BIOMARKERS

    if not (base / "bin/testRun").exists():
        raise FileNotFoundError(
            f"{base / 'bin/testRun'} not found. The ABM is invoked relative "
            f"to its working directory ({base.resolve()}). Run from the "
            f"directory holding bin/, or use prepare_workdir to set one up."
        )

    if biomarkers.exists():
        biomarkers.unlink()

    env = os.environ.copy()
    if device is not None:
        env["CUDA_VISIBLE_DEVICES"] = str(device)

    # The config path is resolved before we change directory, so callers
    # can pass a path relative to their own cwd
    config_arg = Path(config_file)
    if workdir is not None and not config_arg.is_absolute():
        config_arg = config_arg.resolve()



    with open(stdout_file_name, "w") as stdout_file, open(stderr_file_name, "w") as stderr_file:
        status = subprocess.call([
            "./bin/testRun",
            "--numticks",
            str(numticks),
            "--config",
            str(config_arg),
            "--wxw",
            "0.6",
            "--wyw",
            "0.6",
            "--wzw",
            "0.6"
        ], stdout=stdout_file, stderr=stderr_file, cwd=str(base), env=env)

    if status != 0:
        tail = ""
        if stderr_file_name.exists():
            lines = stderr_file_name.read_text(errors="replace").splitlines()
            tail = "\n".join(lines[-15:])
        raise RuntimeError(
            f"./bin/testRun exited with status {status} for config "
            f"{config_arg} ({numticks} ticks, workdir {base}).\n"
            f"Last lines of {stderr_file_name}:\n{tail}"
        )

    if not biomarkers.exists():
        raise RuntimeError(
            f"./bin/testRun reported success but wrote no {biomarkers}. "
            f"Check {stdout_file_name} and {stderr_file_name}."
        )

    return biomarkers
