import os
import re
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
    parameters.xlsx, so the two can be cross-validated.
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
    "Parameter Name" column of parameters.xlsx.
    """
    lines = load_template_lines(template_file)
    return [path for _, path in find_variable_line_entries(lines)]

def create_sample_file(
    parameter_list: list[float],
    template_file: Path,
    out_path: Path,
    parameter_names: list[str] | None = None,
    ) -> None:
    """
    Create a sample file (JSON template) based on the given parameters.
    All other parts of the template (world_init, chemistry) are copied through
    unchanged
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
    
    if out_path.exists():
        os.remove(out_path)
    
    with open(out_path, "w") as f:
        f.writelines(new_lines)

def extract_output_metrics(output_file: Path) -> dict[str, float]:
    """
    Extract output metrics from the ABM output_biomarkers.csv file.
    Modify days and biomarkers here as needed (whatever you want to optimize against)
    """
    df = pd.read_csv(output_file)
    day3_tick = 132
    day6_tick = 264
    day3_row = df[df["clock"] == day3_tick]
    day6_row = df[df["clock"] == day6_tick]
    return {
        "day_3_collagen": day3_row["Collagen (ug)"].values[0],
        "day_3_cells": day3_row["Total Cells"].values[0],
        "day_6_collagen": day6_row["Collagen (ug)"].values[0],
        "day_6_cells": day6_row["Total Cells"].values[0],
    }

def run_ABM(config_file: Path) -> None:
    """
    Run the ABM with some given configuration file.
    """
    stdout_file_name = "output/stdout.txt"
    stderr_file_name = "output/stderr.txt"
    with open(stdout_file_name, "w") as stdout_file, open(stderr_file_name, "w") as stderr_file:  
        subprocess.call([
            "./bin/testRun",
            "--numticks",
            "289", 
            "--inputfile",
            str(config_file),
            "--wxw",
            "0.6",
            "--wyw",
            "0.6",
            "--wzw",
            "0.6"
        ], stdout=stdout_file, stderr=stderr_file)