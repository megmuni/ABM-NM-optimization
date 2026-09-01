from scipy.optimize import minimize
from typing import Any
from pathlib import Path
import shutil
import os
import numpy as np
import pandas as pd
import subprocess
import csv
import argparse 

from ABM import create_sample_file, run_ABM

# ==========================
# CONSTANTS
# TODO: Move these into a config file.
SNAPSHOT_INTERVAL = 5
TICKS_PER_DAY = 48

CONFIG_TEMPLATE = Path("configFiles/simulation_config.template.json")

CONDITION_TO_GROUP = {
    "high": "config_scaffold_High",
    "low": "config_scaffold_Low",
}
# ==========================

# ==========================
# Command line arguments
# n: number of parameters to optimize
# method: method to use for parameter importance ranking, must
# match the column name in the excel file
# "parameters.xlsx"
# ==========================
parser = argparse.ArgumentParser(description='Run ABM optimization.')

parser.add_argument('--n', type=int, default=5, help='Number of parameters to optimize')
#parser.add_argument('--method', type=str, default='Random Forest', help='Method to use for parameter importance ranking')
parser.add_argument('--snapshots', type=bool, default=False, help='Save snapshots of the biomarker values')

args = parser.parse_args() 
n = args.n 
#method = args.method
# ==========================
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

def extract_varying_param_indices() -> list:
    """
    Extract the indices of every parameter that should be optimized --
    i.e. every row in parameters.xlsx whose "Vary?" column is
    blank/NaN, matching the same convention used in the upstream
    sampling/RF pipeline (blank = vary, any non-blank value e.g. "N" =
    don't vary, use the default).
    """
    df = pd.read_excel("parameters.xlsx")
 
    # Same robust "blank vs non-blank" check used upstream: treat any
    # non-empty value in "Vary?" as "don't vary", rather than matching the
    # literal string "N" (robust to whitespace/case/type quirks in Excel).
    vary_mask = df["Vary?"].apply(lambda v: not (pd.notna(v) and str(v).strip() != ""))
 
    varying_rows = df[vary_mask]
 
    if varying_rows.empty:
        raise ValueError('No parameters found with a blank "Vary?" column -- nothing to optimize.')
 
    param_nums = [int(p) for p in varying_rows["Parameter Number"].tolist()]
 
    print(f"Optimizing {len(param_nums)} parameters (Vary? blank): {param_nums}")
 
    return param_nums

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

def extract_biomarkers_at_tick(csv_path: str, tick: int) -> dict[str, float]:
    """
    Read Output_Biomarkers.csv and return the biomarker values at the row
    whose clock column equals `tick`
    """
    df = pd.read_csv(csv_path)
    row = df[df[CLOCK_COL] == tick]
    if row.empty:
        raise ValueError(f"No row found in {csv_path} with {CLOCK_COL} == {tick}")
 
    return {
        "aggrecan": float(row[AGGRECAN_COL].values[0]),
        "total_cells": float(row[TOTAL_CELLS_COL].values[0]),
        "cell_viability": float(row[CELL_VIABILITY_COL].values[0]),
        "percent_diff": float(row[PERCENT_DIFF_COL].values[0]),
    }

# Column names taken directly from Output_Biomarkers.csv header row
CLOCK_COL = "clock (30 min)"
AGGRECAN_COL = "Aggrecan (ug)"
TOTAL_CELLS_COL = "Total Cells"
CELL_VIABILITY_COL = "Viability Rate(%)"
PERCENT_DIFF_COL = "Differentiation (%)"

def save_snapshot(nfeval: int, config: str):
    """
    Save snapshot for this configuration to a CSV file.
    """
    snapshot_data_day_7 = extract_row_as_dict('output/Output_Biomarkers.csv', TICKS_PER_DAY * 7)
    snapshot_data_day_21 = extract_row_as_dict('output/Output_Biomarkers.csv', TICKS_PER_DAY * 21)
    snapshot_data_day_7_with_nfeval = dict(snapshot_data_day_7)
    snapshot_data_day_7_with_nfeval['nfeval'] = nfeval
    snapshot_data_day_7_with_nfeval['config'] = config

    snapshot_data_day_21_with_nfeval = dict(snapshot_data_day_21)
    snapshot_data_day_21_with_nfeval['nfeval'] = nfeval
    snapshot_data_day_21_with_nfeval['config'] = config
    fieldnames = list(snapshot_data_day_7.keys()) + ['nfeval'] + ['config']
    with open('output/snapshots/day_snapshots.csv', 'a', newline='') as snapshots_file:
        writer = csv.DictWriter(snapshots_file, fieldnames=fieldnames)
        if snapshots_file.tell() == 0:
            writer.writeheader()
        writer.writerow(snapshot_data_day_7_with_nfeval)
        writer.writerow(snapshot_data_day_21_with_nfeval)

def run_with_scaffold(condition: str, Y: np.ndarray, experimental_df: pd.DataFrame, config_index: int, param_values: list, 
                      param_names: list, num_iters: int = 3):
    """
    Run the ABM simulation with the specified scaffold configuration 
    # ("high" or "low")
    """
    # Collect errors for each run
    # modify this list based on what experimental data you have available,
    # at which timepoints
    cellviability_day7_errors = []
    percentdiff_day7_errors = []
    cellviability_day21_errors = []
    aggrecan_day21_errors = []
    percentdiff_day21_errors = []
    
    sample_config_path = Path("configFiles/simulation_config_sample.json")

    for iter in range(num_iters):
        print(f"Running iteration {iter + 1} for config {config_file}")
        with open(stdout_file_name, 'a') as stdout_file:
            with open(stderr_file_name, 'a') as stderr_file:
                stdout_file.write("\n\n******************************\n*** MODEL EXECUTION #" + str(Nfeval) + " ***\n******************************\n")
                stderr_file.write("\n\n******************************\n*** MODEL EXECUTION #" + str(Nfeval) + " ***\n******************************\n")
                create_sample_file(param_values, CONFIG_TEMPLATE, sample_config_path,
                            parameter_names=param_names, condition=condition)
                run_ABM(sample_config_path)

        # After each run, read output and calculate error
        with open('output/Output_Biomarkers.csv', 'rt') as f:
            temp = csv.reader(f)
            temp = list(temp)

            group_name = Path(config_file).stem
            day7_cellviability = experimental_df.loc[(group_name, 168), 'small_scaffold_cell_viability']
            day21_cellviability = experimental_df.loc[(group_name, 504), 'small_scaffold_cell_viability']
            #day7_aggrecan= experimental_df.loc[(group_name, 168), 'small_scaffold_aggrecan_pg']
            day21_aggrecan = experimental_df.loc[(group_name, 504), 'small_scaffold_aggrecan_ug']
            day7_percentdiff = experimental_df.loc[(group_name, 168), 'small_scaffold_percent_diff']
            day21_percentdiff = experimental_df.loc[(group_name, 504), 'small_scaffold_percent_diff']

            print(f"Day 7: aggrecan={temp[TICKS_PER_DAY * 7][7} total cells={temp[TICKS_PER_DAY * 7][8]} cell viability={temp[TICKS_PER_DAY * 7][20]} % diff={temp[TICKS_PER_DAY * 7][21]}")
            print(f"Day 21: aggrecan={temp[TICKS_PER_DAY * 21][7]} total cells={temp[TICKS_PER_DAY * 21][8]} cell viability={temp[TICKS_PER_DAY * 21][20]} % diff={temp[TICKS_PER_DAY * 21][21]}")

            cellviability_day7_errors.append(error(day7_cellviability, float(temp[TICKS_PER_DAY * 7][20])))
            percentdiff_day7_errors.append(error(day7_percentdiff, float(temp[TICKS_PER_DAY * 7][21])))
            
            aggrecan_day21_errors.append(error(day21_aggrecan, float(temp[TICKS_PER_DAY * 21][7])))
            cellviability_day21_errors.append(error(day21_cellviability, float(temp[TICKS_PER_DAY * 21][20])))
            percentdiff_day21_errors.append(error(day21_percentdiff, float(temp[TICKS_PER_DAY * 21][21])))

    # Assign the mean error over all runs
    base = config_index * 5
    Y[base + 0] = np.mean(cellviability_day7_errors)
    Y[base + 1] = np.mean(percentdiff_day7_errors)
    Y[base + 2] = np.mean(aggrecan_day21_errors)
    Y[base + 3] = np.mean(cellviability_day21_errors)
    Y[base + 4] = np.mean(percentdiff_day21_errors)

def ABM(x):

    global Nfeval
    global experimental_df_indexed
    # Put sampled parameters into text files
    sam = np.asarray(temp_sample)
    for idx, param in enumerate(params):
        sam[param] = x[idx]
    param_values = sam.tolist()

    Y = np.zeros(10)

    run_with_scaffold("high", Y, experimental_df_indexed, config_index=0, 
                      param_values=param_values, param_names=param_names)

    if args.snapshots:
        save_snapshot(Nfeval, "high")
        if Nfeval % SNAPSHOT_INTERVAL == 0:
            shutil.copy('output/Output_Biomarkers.csv', f'output/snapshots/biomarker_csvs/snapshots_{Nfeval}_high.csv')

    run_with_scaffold("low", Y, experimental_df_indexed, config_index=1,
                       param_values=param_values, param_names=param_names)

    if args.snapshots:
        save_snapshot(Nfeval, "low")
        if Nfeval % SNAPSHOT_INTERVAL == 0:
            shutil.copy('output/Output_Biomarkers.csv', f'output/snapshots/biomarker_csvs/snapshots_{Nfeval}_low.csv')

    # Dynamically create string based on the number of parameters
    format_str = formatted_string(Nfeval, x, Y)     
    print(format_str)
    Nfeval += 1 # Increment function evaluation count
    print(Y)

    return np.sum(Y) #SSE

if __name__ == "__main__":
    # Create parameter names
    df = pd.read_excel("parameters.xlsx")
    numpar = len(df)
    param_names = df["Parameter Name"].tolist()

    # Update array with selected parameters
    #params = extract_n_params(method=method, n=n)
    #params = [i for i in range(numpar)]
    #params = [7, 11, 15, 20, 21, 30, 35, 36, 37, 38] #testing with 10 important params
    # Parameters to optimize: every row with a blank "Vary?" column.
    params = extract_varying_param_indices(df)
    
    # Read bounds
    bounds = df[["Lower", "Upper"]].to_numpy()

    # Read default values
    temp_sample_1 = df[["Mean"]].to_numpy()
    global temp_sample
    temp_sample = np.reshape(temp_sample_1,numpar)

    # Choose specific parameters
    names_s = list( names[i] for i in params )
    print(names_s)
    bounds_s = bounds[np.array(params)]
    print(bounds_s)
    default_s = temp_sample[np.array(params)]
    print(default_s)

    # Open files
    stdout_file_name = "output/SensitivityAnalysis/stdout.txt"
    stderr_file_name = "output/SensitivityAnalysis/stderr.txt"
    open(stdout_file_name, 'w').close()
    open(stderr_file_name, 'w').close()

    # Create snapshots file
    if args.snapshots:
        os.makedirs('output/snapshots', exist_ok=True)
        open('output/snapshots/day_snapshots.csv', 'w').close()
        os.makedirs('output/snapshots/biomarker_csvs', exist_ok=True)

    # # # Run optimization # # #
    Nfeval = 1

    # Construct an initial simplex
    selected_init = construct_simplex(bounds, params)

    # Generate experimental values
    global experimental_df_indexed
    experimental_df = extract_small_scaffold_experimental(Path("experimental_config.csv"))
    experimental_df_indexed = experimental_df.set_index(['group', 'time_hour'])

    print("Experimental data:")
    print(experimental_df_indexed)

    result = minimize(
        ABM, 
        default_s, 
        method='nelder-mead', 
        tol = 1e-4, 
        bounds = bounds_s, 
        options={
            'maxiter': 50, 
            'disp': True, 
            'initial_simplex': selected_init, 
            'return_all': True
        }
    )

    result.x, result.fun

    print(result.x)
    print(result.fun)
    print(result.message)

