from scipy.optimize import minimize
import numpy as np
import pandas as pd
import subprocess
import csv
import os

def ABM(x):

    global Nfeval

    # Put sampled parameters into text files
    sam = np.asarray(temp_sample)
    sam[72] = x[0]
    sam[40] = x[1]
    sam[47] = x[2] 
    sam[21] = x[3]
    sam[33] = x[4]

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


        print('{0:4d}   {1: 3.6f}  {2: 3.6f}  {3: 3.6f} {4: 3.6f}  {5: 3.6f}    {6: 3.6f}'.format(Nfeval, x[0], x[1], x[2], x[3], x[4], np.sum(Y)))
        Nfeval += 1

    return np.sum(Y) #SSE

#####

# Create parameter names
numpar = 75 # total number of parameters
names = ["" for j in range(numpar)]

p1 = 72
p2 = 40
p3 = 47
p4 = 21
p5 = 33

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
par_s = 5 # of parameters chosen by sensitivity analysis
names_s = list( names[i] for i in [p1,p2,p3,p4,p5] )
print(names_s)
bounds_s = bounds[np.array([p1,p2,p3,p4,p5])]
print(bounds_s)
default_s = temp_sample[np.array([p1,p2,p3,p4,p5])]
print(default_s)

# Open files
stdout_file_name = "output/SensitivityAnalysis/stdout.txt"
stderr_file_name = "output/SensitivityAnalysis/stderr.txt"
open(stdout_file_name, 'w').close()
open(stderr_file_name, 'w').close()

# # # Run optimization # # #
Nfeval = 1

init = [[2,0,0,2.5,1],[200, 0, 0,2.5,1],[2, 50, 0,2.5,1],[2, 0, 9630,2.5,1],[2,0,0,100,1],[200,50,9639,100,100]] # initial values

result = minimize(ABM, default_s, method='nelder-mead', tol = 1e-4, bounds = bounds_s, options={'maxiter': 50, 'disp': True, 'initial_simplex': init, 'return_all': True}) 

result.x, result.fun

print(result.x)
print(result.fun)
print(result.message)
