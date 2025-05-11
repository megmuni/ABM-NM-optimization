from scipy.optimize import minimize
import numpy as np
import pandas as pd
import subprocess
import csv
import argparse 

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
parser.add_argument('--test', type=bool, default=False, help='Test mode (True/False)')

args = parser.parse_args() 
n = args.n 
method = args.method
# ==========================
# ==========================

def construct_simplex(bounds: np.ndarray, selected_params: list) -> np.ndarray:
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
    
    # convert to numpy array
    return np.array(init_simplex)

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

def test_ABM(x):
    """
    Unused, but kept for testing purposes.
    The ABM function is very resource intensive, so this
    function allows for a quick test for size and shape matching.
    Test function to check if the parameters are being passed correctly.
    """
    
    global Nfeval
    sam = np.asarray(temp_sample)
    for idx, param in enumerate(params):
        sam[param] = x[idx]

    format_str = formatted_string(Nfeval, x, np.zeros((6,4)))     
    print(format_str)
        
    return 0

def ABM(x):

    global Nfeval

    # Put sampled parameters into text files
    sam = np.asarray(temp_sample)
    for idx, param in enumerate(params):
        sam[param] = x[idx]

    # Put sampled parameters into text files
    np.savetxt("Sample.txt", [sam], delimiter='\t')

    Y = np.zeros((6,4))

    for i in range(3):

        # Run model
        with open(stdout_file_name, 'a') as stdout_file:
            with open(stderr_file_name, 'a') as stderr_file:
                stdout_file.write("\n\n******************************\n*** MODEL EXECUTION #" + str(Nfeval) + " ***\n******************************\n")
                stderr_file.write("\n\n******************************\n*** MODEL EXECUTION #" + str(Nfeval) + " ***\n******************************\n")
                subprocess.call(["./bin/testRun", "--numticks", "289" , "--inputfile" , "configFiles/config_Scaffold_GH2.txt", "--wxw", "0.6", "--wyw", "0.6", "--wzw", "0.6"], stdout = stdout_file, stderr = stderr_file)

        # Save output
        with open('output/Output_Biomarkers.csv', 'rt') as f:
            temp = csv.reader(f)
            temp = list(temp)

            # # Day 3
            Y[0][i] = ((float(temp[144][16]) + float(temp[144][17]) - 90)/max(float(temp[144][16]) + float(temp[144][17]),90))**2 # Fibroblasts
            Y[1][i] = ((float(temp[144][8]) - 64736.8)/max(float(temp[144][8]),64736.8))**2 # Collagen
            
            # # Day 6
            Y[0][i] = ((float(temp[288][16]) + float(temp[288][17]) - 90)/max(float(temp[288][16]) + float(temp[288][17]),90))**2 # Fibroblasts
            Y[1][i] = ((float(temp[288][8]) - 42785)/max(float(temp[288][8]),42785))**2 # Collagen

            # Validation
            # Y[0][i] = ((float(temp[4][18]) + float(temp[4][21]) - 3981)/max(float(temp[4][18]) + float(temp[4][21]),3981))**2 # Fibroblasts
            # Y[1][i] = ((float(temp[4][9]) + float(temp[4][10]) + float(temp[4][11]) - 80860)/max(float(temp[4][9]) + float(temp[4][10]) + float(temp[4][11]),80860))**2 # Collagen

        # Run model
        with open(stdout_file_name, 'a') as stdout_file:
            with open(stderr_file_name, 'a') as stderr_file:
                stdout_file.write("\n\n******************************\n*** MODEL EXECUTION #" + str(Nfeval) + " ***\n******************************\n")
                stderr_file.write("\n\n******************************\n*** MODEL EXECUTION #" + str(Nfeval) + " ***\n******************************\n")
                subprocess.call(["./bin/testRun", "--numticks", "289" , "--inputfile" , "configFiles/config_Scaffold_GH5.txt", "--wxw", "0.6", "--wyw", "0.6", "--wzw", "0.6"], stdout = stdout_file, stderr = stderr_file)

        # Save output
        with open('output/Output_Biomarkers.csv', 'rt') as f:
            temp = csv.reader(f)
            temp = list(temp)

            # # Day 3
            Y[2][i] = ((float(temp[144][16]) + float(temp[144][17]) - 90)/max(float(temp[144][16]) + float(temp[144][17]),90))**2 # Fibroblasts
            Y[3][i] = ((float(temp[144][8]) - 64736.8)/max(float(temp[144][8]),64736.8))**2 # Collagen
            
            # # Day 6
            Y[2][i] = ((float(temp[288][16]) + float(temp[288][17]) - 90)/max(float(temp[288][16]) + float(temp[288][17]),90))**2 # Fibroblasts
            Y[3][i] = ((float(temp[288][8]) - 42785)/max(float(temp[288][8]),42785))**2 # Collagen

            # Validation
            # Y[0][i] = ((float(temp[4][18]) + float(temp[4][21]) - 3981)/max(float(temp[4][18]) + float(temp[4][21]),3981))**2 # Fibroblasts
            # Y[1][i] = ((float(temp[4][9]) + float(temp[4][10]) + float(temp[4][11]) - 80860)/max(float(temp[4][9]) + float(temp[4][10]) + float(temp[4][11]),80860))**2 # Collagen

        # Run model
        with open(stdout_file_name, 'a') as stdout_file:
            with open(stderr_file_name, 'a') as stderr_file:
                stdout_file.write("\n\n******************************\n*** MODEL EXECUTION #" + str(Nfeval) + " ***\n******************************\n")
                stderr_file.write("\n\n******************************\n*** MODEL EXECUTION #" + str(Nfeval) + " ***\n******************************\n")
                subprocess.call(["./bin/testRun", "--numticks", "289" , "--inputfile" , "configFiles/config_Scaffold_GH10.txt", "--wxw", "0.6", "--wyw", "0.6", "--wzw", "0.6"], stdout = stdout_file, stderr = stderr_file)

        # Save output
        with open('output/Output_Biomarkers.csv', 'rt') as f:
            temp = csv.reader(f)
            temp = list(temp)

            # # Day 3
            Y[4][i] = ((float(temp[144][16]) + float(temp[144][17]) - 90)/max(float(temp[144][16]) + float(temp[144][17]),90))**2 # Fibroblasts
            Y[5][i] = ((float(temp[144][8]) - 64736.8)/max(float(temp[144][8]),64736.8))**2 # Collagen
            
            # # Day 6
            Y[4][i] = ((float(temp[288][16]) + float(temp[288][17]) - 90)/max(float(temp[288][16]) + float(temp[288][17]),90))**2 # Fibroblasts
            Y[5][i] = ((float(temp[288][8]) - 42785)/max(float(temp[288][8]),42785))**2 # Collagen

            # Validation
            # Y[0][i] = ((float(temp[4][18]) + float(temp[4][21]) - 3981)/max(float(temp[4][18]) + float(temp[4][21]),3981))**2 # Fibroblasts
            # Y[1][i] = ((float(temp[4][9]) + float(temp[4][10]) + float(temp[4][11]) - 80860)/max(float(temp[4][9]) + float(temp[4][10]) + float(temp[4][11]),80860))**2 # Collagen

        # Dynamically create string based on the number of parameters
        format_str = formatted_string(Nfeval, x, Y)     
        #print('{0:4d}   {1: 3.6f}  {2: 3.6f}  {3: 3.6f} {4: 3.6f}  {5: 3.6f}    {6: 3.6f}'.format(Nfeval, x[0], x[1], x[2], x[3], x[4], np.sum(Y)))
        print(format_str)
        NFeval += 1

    return np.sum(Y) #SSE

# =================================

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

# # # Run optimization # # #
Nfeval = 1

# Construct an initial simplex
selected_init = construct_simplex(bounds, params)


result = minimize(
    ABM if not args.test else test_ABM, 
    default_s, 
    method='nelder-mead', 
    tol = 1e-4, 
    bounds = bounds_s, 
    options={
        'maxiter': 50 if not args.test else 5, 
        'disp': True, 
        'initial_simplex': selected_init, 
        'return_all': True
    }
)

result.x, result.fun

print(result.x)
print(result.fun)
print(result.message)

