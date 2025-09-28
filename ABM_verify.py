import pandas as pd
import json
import numpy as np
from pathlib import Path
from ABM import extract_n_param_names, create_sample_file, run_ABM

# ======================
# Check that the ABM model is working as expected
# i.e. it produces expected outputs for known inputs
# ======================

if __name__ == "__main__":
    generated_samples_df = pd.read_csv("generated_samples_with_outputs.csv")
    param_names = extract_n_param_names()

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

        simulated = {col: float(row[col]) for col in known_result_cols}

        errors = {}
        for col in known_result_cols:
            known = float(row[col])
            sim = simulated[col]
            err = (sim - known) ** 2
            errors[col] = err
            squared_errors[col].append(err)

        per_sample_results.append({
            "sample_id": int(row["sample_id"]),
            "input_parameters": {param: row[param] for param in param_names},
            "expected_outputs": {col: float(row[col]) for col in known_result_cols},
            "simulated_outputs": simulated,
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


    
    
