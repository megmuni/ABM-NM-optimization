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

# ==========================
# CONSTANTS
# TODO: Move these into a config file.
FIBROBLASTS = 90
DAY_3_COLLAGEN = 64736.8
DAY_6_COLLAGEN = 42785
SNAPSHOT_INTERVAL = 5
TICKS_PER_DAY = 44
# ==========================

# ==========================
# Command line arguments
# n: number of parameters to optimize
# method: method to use for parameter importance ranking, must
# match the column name in the excel file
# "Sensitivity Analysis.xlsx"
# ==========================
parser = argparse.ArgumentParser(description='Run ABM optimization.')

parser.add_argument('--n', type=int, default=5, help='Number of parameters to optimize')
parser.add_argument('--method', type=str, default='Random Forest', help='Method to use for parameter importance ranking')
parser.add_argument('--snapshots', type=bool, default=False, help='Save snapshots of the biomarker values')

args = parser.parse_args() 
n = args.n 
method = args.method
# ==========================
# ==========================


def small_scaffold_adjustment_cells(value: float) -> float:
    """
    Adjusts the value for the small scaffold experimental data.
    The adjustment is based on the formula:
    ((0.6)^3 / 1000 / 0.3 ) * value
    """
    return ((0.6 ** 3) / 1000 / 0.3) * value

def small_scaffold_adjustment_collagen(value: float) -> float:
    """
    Adjusts the value for the small scaffold experimental data.
    The adjustment is based on the formula:
    ((0.6)^3 / 300 ) * value * 10^6
    """
    return ((0.6 ** 3) / 300) * value * 10 ** 6

def extract_small_scaffold_experimental(file_path: Path) -> pd.DataFrame:
    """
    Extracts experimental data from CSV file corresponding to the small scaffold.
    """
    df = pd.read_csv(file_path)
    df.columns = [c.strip() for c in df.columns]

    # Calculate average for each config (GH2, GH5, GH10)
    df["picogreen_cells"] = df["picogreen_cells"].astype(float)
    df["sircol_collagen_ug"] = df["sircol_collagen_ug"].astype(float)
    df["bradford_protein_ug_per_ml"] = df["bradford_protein_ug_per_ml"].astype(float)

    # Add a column for live cells
    df["live_cells"] = df["picogreen_cells"] * (df["cell_viability"])

    # Calculate the average values for each group and time point
    averages = df.groupby(["group", "time_hour"]).mean().reset_index()

    # Apply scaffold adjustment to the mean live_cells
    averages["small_scaffold_cell_avg"] = averages["live_cells"].apply(small_scaffold_adjustment_cells)

    # Collagen adjustment (convert from ug to pg)
    averages["small_scaffold_collagen_pg"] = averages["sircol_collagen_ug"].apply(
        lambda x: small_scaffold_adjustment_collagen(x)
    )
    
    return averages

def error(expected, actual):
    """
    Calculate SSE between expected and actual values.
    """
    return ((expected - actual) / max(expected, actual)) ** 2

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

def extract_n_params(method="Random Forest", n=5) -> list:
    """
    Extract the n most important parameters from the given method.
    Default: Random Forest, n = 5
    """
    # read senstivity analysis data and drop unranked parameters
    df = pd.read_excel("Sensitivity Analysis.xlsx")
    print(df)
    df_filtered = df.dropna(subset=[method])
    
    # sort by the column corresponding to the method
    df_sorted = df_filtered.sort_values(by=method, ascending=True)
    
    # get the top n parameters
    top_n_params = df_sorted.head(n)
    
    if len(top_n_params) < n:
        raise ValueError(f"Not enough parameters found for method: {method}. Found: {top_n_params.shape[0]}, expected: {n}.")

    # extract the parameter numbers as a list
    param_nums = [int(param) for param in top_n_params["Parameter Number"].tolist()]
    
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

def save_snapshot(nfeval: int, config: str):
    """
    Save snapshot for this configuration to a CSV file.
    """
    snapshot_data_day_3 = extract_row_as_dict('output/Output_Biomarkers.csv', TICKS_PER_DAY * 3)
    snapshot_data_day_6 = extract_row_as_dict('output/Output_Biomarkers.csv', TICKS_PER_DAY * 6)
    snapshot_data_day_3_with_nfeval = dict(snapshot_data_day_3)
    snapshot_data_day_3_with_nfeval['nfeval'] = nfeval
    snapshot_data_day_3_with_nfeval['config'] = config

    snapshot_data_day_6_with_nfeval = dict(snapshot_data_day_6)
    snapshot_data_day_6_with_nfeval['nfeval'] = nfeval
    snapshot_data_day_6_with_nfeval['config'] = config
    fieldnames = list(snapshot_data_day_3.keys()) + ['nfeval'] + ['config']
    with open('output/snapshots/day_snapshots.csv', 'a', newline='') as snapshots_file:
        writer = csv.DictWriter(snapshots_file, fieldnames=fieldnames)
        if snapshots_file.tell() == 0:
            writer.writeheader()
        writer.writerow(snapshot_data_day_3_with_nfeval)
        writer.writerow(snapshot_data_day_6_with_nfeval)

def run_with_scaffold(config_file: str, Y: np.ndarray, experimental_df: pd.DataFrame, i: int, num_iters: int = 3):
    """
    Run the ABM simulation with the specified scaffold configuration file.
    """
    # Collect errors for each run
    fibroblast_day3_errors = []
    collagen_day3_errors = []
    fibroblast_day6_errors = []
    collagen_day6_errors = []

    for iter in range(num_iters):
        print(f"Running iteration {iter + 1} for config {config_file}")
        with open(stdout_file_name, 'a') as stdout_file:
            with open(stderr_file_name, 'a') as stderr_file:
                stdout_file.write("\n\n******************************\n*** MODEL EXECUTION #" + str(Nfeval) + " ***\n******************************\n")
                stderr_file.write("\n\n******************************\n*** MODEL EXECUTION #" + str(Nfeval) + " ***\n******************************\n")
                subprocess.call(["./bin/testRun", "--numticks", "289" , "--inputfile" , config_file, "--wxw", "0.6", "--wyw", "0.6", "--wzw", "0.6"], stdout = stdout_file, stderr = stderr_file)

        # After each run, read output and calculate error
        with open('output/Output_Biomarkers.csv', 'rt') as f:
            temp = csv.reader(f)
            temp = list(temp)

            group_name = Path(config_file).stem
            day3_fibroblasts = experimental_df.loc[(group_name, 72), 'small_scaffold_cell_avg']
            day6_fibroblasts = experimental_df.loc[(group_name, 144), 'small_scaffold_cell_avg']
            day3_collagen = experimental_df.loc[(group_name, 72), 'small_scaffold_collagen_pg']
            day6_collagen = experimental_df.loc[(group_name, 144), 'small_scaffold_collagen_pg']

            print(f"Day 3: collagen={temp[TICKS_PER_DAY * 3][8]} activated={temp[TICKS_PER_DAY * 3][16]} fibroblasts={temp[TICKS_PER_DAY * 3][17]}")
            print(f"Day 6: collagen={temp[TICKS_PER_DAY * 6][8]} activated={temp[TICKS_PER_DAY * 6][16]} fibroblasts={temp[TICKS_PER_DAY * 6][17]}")

            fibroblast_day3_errors.append(error(day3_fibroblasts, float(temp[TICKS_PER_DAY * 3][16]) + float(temp[TICKS_PER_DAY * 3][17])))
            collagen_day3_errors.append(error(day3_collagen, float(temp[TICKS_PER_DAY * 3][8])))
            fibroblast_day6_errors.append(error(day6_fibroblasts, float(temp[TICKS_PER_DAY * 6][16]) + float(temp[TICKS_PER_DAY * 6][17])))
            collagen_day6_errors.append(error(day6_collagen, float(temp[TICKS_PER_DAY * 6][8])))

    # Assign the mean error over all runs
    Y[0][i] = np.mean(fibroblast_day3_errors)
    Y[1][i] = np.mean(collagen_day3_errors)
    Y[2][i] = np.mean(fibroblast_day6_errors)
    Y[3][i] = np.mean(collagen_day6_errors)

def ABM(x):

    global Nfeval

    global experimental_df_indexed
    # Put sampled parameters into text files
    sam = np.asarray(temp_sample)
    for idx, param in enumerate(params):
        sam[param] = x[idx]

    # Put sampled parameters into text files
    np.savetxt("Sample.txt", [sam], delimiter='\t')

    Y = np.zeros((12, 4))

    for i in range(3):
        
        run_with_scaffold("configFiles/config_Scaffold_GH2.txt", Y, experimental_df_indexed, i)

        if args.snapshots:
            save_snapshot(Nfeval, "config_Scaffold_GH2")
            if Nfeval % SNAPSHOT_INTERVAL == 0:
                shutil.copy('output/Output_Biomarkers.csv', f'output/snapshots/biomarker_csvs/snapshots_{Nfeval}_config_Scaffold_GH2.csv')

        run_with_scaffold("configFiles/config_Scaffold_GH5.txt", Y, experimental_df_indexed, i)

        if args.snapshots:
            save_snapshot(Nfeval, "config_Scaffold_GH5")
            if Nfeval % SNAPSHOT_INTERVAL == 0:
                shutil.copy('output/Output_Biomarkers.csv', f'output/snapshots/biomarker_csvs/snapshots_{Nfeval}_config_Scaffold_GH5.csv')

        run_with_scaffold("configFiles/config_Scaffold_GH10.txt", Y, experimental_df_indexed, i)

        if args.snapshots:
            save_snapshot(Nfeval, "configFiles/config_Scaffold_GH10")
            if Nfeval % SNAPSHOT_INTERVAL == 0:
                shutil.copy('output/Output_Biomarkers.csv', f'output/snapshots/biomarker_csvs/snapshots_{Nfeval}_config_Scaffold_GH10.csv')

        # Dynamically create string based on the number of parameters
        format_str = formatted_string(Nfeval, x, Y)     
        print(format_str)
        Nfeval += 1

    return np.sum(Y) #SSE

if __name__ == "__main__":
    # Create parameter names
    numpar = 75 # total number of parameters
    names = ["" for j in range(numpar)]

    # Update array with selected parameters
    params = extract_n_params(method=method, n=n)

    df = pd.read_excel(r'Sensitivity Analysis.xlsx') # read parameter bounds

    for i in range(numpar):
        names[i] = "x" + str(i)

    # Read bounds
    bounds = df[["Lower bound", "Upper bound"]].to_numpy()

    # Read default values
    temp_sample_1 = df[["Default Value"]].to_numpy()
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

