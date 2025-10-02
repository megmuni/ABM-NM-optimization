import pandas as pd
import json
import numpy as np
from pathlib import Path
from ABM import extract_n_param_names, create_sample_file, run_ABM, extract_output_metrics

# ======================
# Check that the ABM model is working as expected
# i.e. it produces expected outputs for known inputs
# ======================

def extract_param_names_from_csv(csv_file: str = "generated_samples_with_outputs.csv") -> list[str]:
    """
    Extract parameter names from the CSV file in the order they appear as columns.
    """
    df = pd.read_csv(csv_file)
    
    metadata_cols = ["sample_id", "param_set_id", "configuration"]
    output_cols = ["day_3_collagen", "day_3_fibroblast", "day_6_collagen", "day_6_fibroblast"]
    
    all_cols = df.columns.tolist()
    
    param_cols = [col for col in all_cols if col not in metadata_cols and col not in output_cols]
    
    return param_cols

if __name__ == "__main__":
    generated_samples_df = pd.read_csv("generated_samples_with_outputs.csv")
    param_names = extract_param_names_from_csv()

    known_result_cols = [
        "day_3_collagen",
        "day_3_fibroblast",
        "day_6_collagen",
        "day_6_fibroblast"
    ]

    per_sample_results = []
    squared_errors = {col: [] for col in known_result_cols}

    for idx, row in generated_samples_df.iterrows():
        param_values = [row[param] for param in param_names]
        create_sample_file(param_values)
        run_ABM(Path("Sample.txt"))

        # Extract simulated results from ABM output
        real_output = extract_output_metrics(Path("output/Output_Biomarkers.csv"))
        errors = {}
        for col in known_result_cols:
            known = float(row[col])
            sim = float(real_output[col])
            err = (sim - known) ** 2
            errors[col] = err
            squared_errors[col].append(err)

        per_sample_results.append({
            "sample_id": int(row["sample_id"]),
            "input_parameters": {param: row[param] for param in param_names},
            "expected_outputs": {col: float(row[col]) for col in known_result_cols},
            "simulated_outputs": {col: float(real_output[col]) for col in known_result_cols},
            "squared_errors": errors
        })

    summary_metrics = {
        "mean_squared_error": {col: float(np.mean(squared_errors[col])) for col in known_result_cols},
        "total_samples": len(per_sample_results)
    }

    results_json = {
        "per_sample_results": per_sample_results,
        "summary_metrics": summary_metrics
    }

    with open("verification_results.json", "w") as f:
        json.dump(results_json, f, indent=2)
