from pathlib import Path
from ABM import create_sample_file, run_ABM, extract_output_metrics
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
input_file = "parameters.xlsx"
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
    Run the full pipeline to generate parameter sets with expected outputs.
    """
    param_df = generate_param_df()
    num_varying = param_df["Vary?"].apply(
    lambda v: not (isinstance(v, str) and v.strip().upper() == "N")).sum()
    rows_list = []
    
    #JSON template containing all of the config info (incl. the bio
    #params to be varied)
    config_template = Path("configFiles/simulation_config.template.json")
    
    param_names = param_df["Parameter Name"].tolist()
    conditions = ["high", "low"]
    
    total_generated = 0
    # set value column to default values
    # param_df["value"] = param_df["Default Value"]
    for i in range(n_samples):
        sample_df = mutate_parameters(param_df)
        param_values = sample_df["value"].tolist()
        
        for condition in conditions:
            sample_config_path = Path("configFiles/simulation_config_sample.json")
            create_sample_file(param_values, config_template, sample_config_path, parameter_names=param_names) #make sample param set to run ABM
        
            run_ABM(sample_config_path) #run the ABM with the sample
            metrics = extract_output_metrics(Path("output/Output_Biomarkers.csv"))
            result_row = {
                "sample_id": i,
                "num_params_varied": num_varying,
                **{name: param_values[j] for j, name in enumerate(param_names)},
                **metrics
            }
            rows_list.append(result_row)
            total_generated += 1

    result_df = pd.DataFrame(rows_list)
    return result_df

if __name__ == "__main__":
    df = generate_samples(n_samples=10)
    df.to_csv("generated_samples_with_outputs.csv", index=False)