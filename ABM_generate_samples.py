from pathlib import Path
from ABM import (
    create_sample_file,
    run_ABM,
    extract_output_metrics,
    CONFIG_TEMPLATE,
    SAMPLE_CONFIG_PATH,
    OUTPUT_BIOMARKERS,
    CONDITIONS,
    PARAMETER_FILE,
    format_output_metrics,
)
import numpy as np
import pandas as pd
from scipy.stats import truncnorm

# ==================================================
# Script to generate samples, including parameter sets
# and expected outputs, for the ABM model.
#
# Sept 2026: Updated to work with the JSON config structure.
# A single simulation_config.template.json (with mutable params
# under "biology") is used to write a merged JSON config file per sample
# instead of the old Sample.txt
# ==================================================

# --- Settings to edit -------------------------------------------------
input_file = PARAMETER_FILE
sheet_name = "Sheet1"
n_par = 67  # Number of parameters
rng_seed = None  # Set an int here for reproducibility, or leave None
# -----------------------------------------------------------------------

rng = np.random.default_rng(rng_seed)

def sample_truncnorm(mean, sigma, lower, upper, size, rng):
    """
    Draw `size` samples from a normal(mean, sigma) distribution truncated
    to [lower, upper], so that bounds [a, b] are a hard min/max.
    """
    a_std = (lower - mean) / sigma
    b_std = (upper - mean) / sigma
    return truncnorm.rvs(a_std, b_std, loc=mean, scale=sigma, size=size, random_state=rng)

def generate_param_df() -> pd.DataFrame:
    """
    Generate a DataFrame of parameters with their bounds and default values from spreadsheet.
    
    Columns:
        'Parameter Number': 1-indexed row/param number
        'Parameter Name': colon-joined path matching the config template JSON
        'Mean': default/mean value; used directly when Vary? is 'N'
        'Lower': lower bound (a) of that param
        'Upper': upper bound (b) of that param
        'SD': std dev, precomputed as (Upper - Lower)/2
        'Type': distribution type: 1=lognormal, 2=discrete lognormal, 3=normal, 4=uniform 
        'Vary?': blank=vary this parameter, 'N'=use Mean for every run
    """
    df = pd.read_excel(input_file, sheet_name=sheet_name)
    df = df.iloc[:n_par].reset_index(drop=True)
    return df

def mutate_parameters(params: pd.DataFrame, mutate_params: list[str]) -> pd.DataFrame:
    """
    Mutate the given parameters in the DataFrame.
    Rows marked 'N' in the 'Vary?' column are held at their Mean value
    instead of randomly sampled
    """
    mutated_params = params.copy()
    for i in range(len(mutated_params)):
        row = mutated_params.iloc[i]
        a = row["Lower"]
        b = row["Upper"]
        sd = row["SD"]
        dist_type = row["Type"]
        mean_val = row["Mean"]
        
        vary_raw = row["Vary?"]
        dont_vary = isinstance(vary_raw, str) and vary_raw.strip().upper() == "N"
        
        if dont_vary:
            mutated_params.loc[i, "value"] = mean_val
            continue
        
        if dist_type == 1:
            # lognormal, truncated to [a, b]
            mu = np.log(mean_val)
            sigma = (np.log(b) - np.log(a)) / 2
            log_r = sample_truncnorm(mu, sigma, np.log(a), np.log(b), 1, rng)
            mutated_params.loc[i, "value"] = np.exp(log_r)[0]
        elif dist_type == 2:
            # discrete lognormal, rounded to nearest 0.5
            mu = np.log(mean_val)
            sigma = (np.log(b) - np.log(a)) / 2
            log_x = sample_truncnorm(mu, sigma, np.log(a), np.log(b), 1, rng)
            x = np.exp(log_x)[0]
            y = np.floor(x)
            frac = x - y
            r = y + 0.5 if frac >= 0.5 else y
            mutated_params.loc[i, "value"] = np.clip(r, a, b)
        elif dist_type == 3:
            # normal, truncated to [a, b]
            r = sample_truncnorm(mean_val, sd, a, b, 1, rng)
            mutated_params.loc[i, "value"] = r[0]
        elif dist_type == 4:
            # uniform between lower and upper bound
            mutated_params.loc[i, "value"] = rng.uniform(a, b)
        else:
            raise ValueError(f"Unknown distribution type '{dist_type}' at parameter {i + 1}")
            
    return mutated_params

def generate_samples(n_samples: int, method: str) -> pd.DataFrame:
    """
    Run the full pipeline to generate parameter sets with expected outputs
    
    Each sample run is once per scaffold condition (high/low MW alginate),
    so n_samples parameter sets produce n_samples * len(CONDITIONS) rows
    
    "Condition" column records which one so that ABM_verify can rebuild the
    exact config that produced each row
    """
    param_df = generate_param_df()
    num_varying = param_df["Vary?"].apply(
    lambda v: not (isinstance(v, str) and v.strip().upper() == "N")).sum()
    rows_list = []
 
    param_names = param_df["Parameter Name"].tolist()
 
    for i in range(n_samples):
        sample_df = mutate_parameters(param_df)
        param_values = sample_df["value"].tolist()
 
        for condition in CONDITIONS:
            # Write the merged JSON config (the sampled biology params plus
            # this sample's scaffold condition)
            create_sample_file(
                param_values,
                CONFIG_TEMPLATE,
                SAMPLE_CONFIG_PATH,
                parameter_names=param_names,
                condition=condition,
            )
 
            run_ABM(SAMPLE_CONFIG_PATH)
            metrics = extract_output_metrics(OUTPUT_BIOMARKERS)
            result_row = {
                "sample_id": i,
                "condition": condition,
                "num_params_varied": num_varying,
                **{name: param_values[j] for j, name in enumerate(param_names)},
                **metrics
            }
            rows_list.append(result_row)
            print(f"Generated sample {i} ({condition}): {metrics}")
 
    result_df = pd.DataFrame(rows_list)
    return result_df
 
if __name__ == "__main__":
    df = generate_samples(n_samples=10)
    df.to_csv("generated_samples_with_outputs.csv", index=False)