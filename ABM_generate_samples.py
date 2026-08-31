from pathlib import Path
from ABM import extract_n_param_names, create_sample_file, run_ABM, extract_output_metrics
import numpy as np
import pandas as pd

# ==================================================
# Script to generate samples, including parameter sets
# and expected outputs, for the ABM model.
#
# Sept 2026: Updated to work with the JSON config structure.
# A single simulation_config.template.json (with mutable params
# under "biology") is used to write a merged JSON config file per sample
# instead of the old Sample.txt
# ==================================================

def mutate_parameters(params: pd.DataFrame, mutate_params: list[str]) -> pd.DataFrame:
    """
    Mutate the given parameters in the DataFrame.
    """
    mutated_params = params.copy()
    for param in mutate_params:
        row = mutated_params[mutated_params["Parameter Name"] == param]
        lower = row["Lower bound"].values[0]
        upper = row["Upper bound"].values[0]
        mutated_value = np.random.uniform(lower, upper)
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
    
    config_template = Path("configFiles/simulation_config.template.json")
    
    total_generated = 0
    # set value column to default values
    param_df["value"] = param_df["Default Value"]
    for i in range(n_samples):
        sample_df = mutate_parameters(param_df, top_params)
        param_values = sample_df["value"].tolist()
        
        sample_config_path = Path("configFiles/simulation_config_sample.json")
        create_sample_file(param_values, config_template, sample_config_path) #make sample param set to run ABM
        
        run_ABM(sample_config_path) #run the ABM with the sample
        metrics = extract_output_metrics(Path("output/Output_Biomarkers.csv"))
        result_row = {
            "sample_id": total_generated,
            "param_set_id": i,
            "num_params_mutated": len(top_params),
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