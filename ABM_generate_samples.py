from pathlib import Path
import subprocess
import random
import os
import numpy as np
import pandas as pd

# ==================================================
# Script to generate samples, including parameter sets
# and expected outputs, for the ABM model.
# ==================================================


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
    # Clear the current sample file
    if os.path.exists("Sample.txt"):
        os.remove("Sample.txt")
    np.savetxt("Sample.txt", parameter_list, delimiter="\t")

def run_ABM(config_file: Path) -> None:
    """
    Run the ABM with some given configuration file.
    """
    stdout_file_name = "output/stdout.txt"
    stderr_file_name = "output/stderr.txt"
    # with open(stdout_file_name, "w") as stdout_file_name, open(stderr_file_name, "w") as stderr_file_name:  
    #     subprocess.call([
    #         "./bin/testRun",
    #         "--numticks",
    #         "289", 
    #         "--inputfile",
    #         str(config_file),
    #         "--wxw",
    #         "0.6",
    #         "--wyw",
    #         "0.6",
    #         "--wzw",
    #         "0.6"
    #     ], stdout=stdout_file_name, stderr=stderr_file_name)

def extract_output_metrics(output_file: Path) -> dict[str, float]:
    """
    Extract output metrics from the ABM output_biomarkers.csv file.
    """
    df = pd.read_csv(output_file)
    day3_tick = 132
    day6_tick = 264
    day3_row = df[df["clock"] == day3_tick]
    day6_row = df[df["clock"] == day6_tick]
    return {
        "day_3_collagen": day3_row["Collagen"].values[0],
        "day_3_fibroblast": day3_row["ActivatedFibroblast"].values[0] + day3_row["Fibroblast"].values[0],
        "day_6_collagen": day6_row["Collagen"].values[0],
        "day_6_fibroblast": day6_row["Fibroblast"].values[0] + day6_row["ActivatedFibroblast"].values[0],
    }

def mutate_parameters(params: pd.DataFrame, mutate_params: list[str]) -> pd.DataFrame:
    """
    Mutate the given parameters in the DataFrame.
    """
    mutated_params = params.copy()
    for param in mutate_params:
        row = mutated_params[mutated_params["Parameter Name"] == param]
        lower = row["Lower bound"].values[0]
        upper = row["Upper bound"].values[0]
        mutated_value = random.uniform(lower, upper)
        mutated_params.loc[mutated_params["Parameter Name"] == param, "value"] = mutated_value
            
    return mutated_params

def generate_param_df() -> pd.DataFrame:
    """
    Generate a DataFrame of parameters with their bounds and default values from spreadsheet.
    """
    df = pd.read_excel("Sensitivity Analysis.xlsx")
    param_df = df[["Class", "Parameter Number", "Parameter Name", "Lower bound", "Upper bound", "Default Value"]]
    return param_df

def get_param_names(df: pd.DataFrame) -> list[str]:
    """
    Get the list of parameter names from the DataFrame.
    """
    return df["Parameter Name"].tolist()

def generate_samples(n_samples: int, method: str) -> pd.DataFrame:
    """
    Run the full pipeline to generate parameter sets with expected outputs.
    """
    param_df = generate_param_df()
    all_params = get_param_names(param_df)
    top_params = extract_n_param_names(method, n=5)
    rows_list = []
    config_files = [
        Path("configFiles/config_Scaffold_GH2.txt"),
        Path("configFiles/config_Scaffold_GH5.txt"),
        Path("configFiles/config_Scaffold_GH10.txt"),
    ]
    total_generated = 0
    # set value column to default values
    param_df["value"] = param_df["Default Value"]
    for i in range(n_samples):
        sample_df = mutate_parameters(param_df, top_params)
        param_values = sample_df["value"].tolist()
        create_sample_file(param_values)
        for config in config_files:
            run_ABM(config)
            metrics = extract_output_metrics(Path("output/Output_Biomarkers.csv"))
            result_row = {
                "sample_id": total_generated,
                "param_set_id": i,
                "configuration": config.stem,
                **{param: sample_df[sample_df["Parameter Name"] == param]["value"].values[0] for param in all_params},
                **metrics
            }
            rows_list.append(result_row)
            total_generated += 1
    result_df = pd.DataFrame(rows_list)
    return result_df

if __name__ == "__main__":
    df = generate_samples(n_samples=10, method="Random Forest")
    df.to_csv("generated_samples_with_outputs.csv", index=False)