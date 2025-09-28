import os
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd

def extract_n_param_names(method="Random Forest", n=5) -> list[str]:
    """
    Extract the top n parameter names based on the given method.
    """
    df = pd.read_excel("Sensitivity Analysis.xlsx")
    df_filtered = df.dropna(subset=[method])
    df_sorted = df_filtered.sort_values(by=method)
    top_n_params = df_sorted.head(n)["Parameter Name"].tolist()
    return top_n_params

def create_sample_file(parameter_list: list[float]) -> None:
    """
    Create a sample file based on the given parameters.
    """
    if os.path.exists("Sample.txt"):
        os.remove("Sample.txt")
    np.savetxt("Sample.txt", parameter_list, delimiter="\t")

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