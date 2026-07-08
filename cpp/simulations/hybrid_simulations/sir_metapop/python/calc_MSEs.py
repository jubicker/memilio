import pandas as pd
import numpy as np
from sklearn.metrics import mean_squared_error
from sklearn.metrics import mean_absolute_percentage_error
import os


def scaled_mean_percentage_error(y_true, y_pred):
    return np.mean(np.abs(y_pred - y_true) / ((np.abs(y_true) + np.abs(y_pred))/2.))


tau_infectious = 7

sample_start = 0
sample_end = 4319

data = {"Sample": [], "I_init": [], "R_init": [], "lambda": [], "max_relation": [], "max_relation_adapted": [], "max_mean": [], "R_0": [], "relation_at_max_mu": [], "max_std": [], "max_std_normalized": [], "MSE_mean": [
], "MSE_var": [], "MAPE_mean": [], "MAPE_var": [], "sMAPE_mean": [], "sMAPE_var": [], "corr_mean": [], "corr_var": []}

base_folder = "/hpc_data/bick_ju/TemporalHybrid/"
smm_folder = base_folder + "SMM/performance_study_SIR/"
moment_folder = base_folder + "Moments/performance_study_SIR/"
save_folder = base_folder + "MSE/SIR/"

for s in range(sample_start, sample_end + 1):
    # Add sample number to data
    data["Sample"].append(s)
    # Read in SMM mean and moments
    SMM_mean = pd.read_csv(smm_folder + "sample_" + str(s) + "/means.csv")
    SMM_moments = pd.read_csv(smm_folder + "sample_" + str(s) + "/moments.csv")
    # Read in ODE mean and moemnts
    ODE_mean = pd.read_csv(moment_folder + "sample_" + str(s) + "/means.csv")
    ODE_moments = pd.read_csv(
        moment_folder + "sample_" + str(s) + "/moments.csv")
    # Read in parameters
    parameters = pd.read_csv(smm_folder + "sample_" +
                             str(s) + "/parameters.csv")
    # Add parameters to data
    data["I_init"].append(parameters["I0"].iloc[0])
    data["R_init"].append(parameters["R0"].iloc[0])
    data["lambda"].append(parameters["lambda"].iloc[0])
    # Add max mean to data
    data["max_mean"].append(SMM_mean["C2"].max())
    # Add max relation to data
    max_relation = np.max(np.sqrt(SMM_moments["M020"])/SMM_mean["C2"])
    data["max_relation"].append(max_relation)
    # Add adapted max relation to data
    max_idx = (np.sqrt(SMM_moments["M020"])/SMM_mean["C2"]).idxmax()
    adapted_max_rel = (0 if pd.isna(max_idx) or (np.sqrt(
        SMM_moments["M020"].iloc[max_idx]) < 1 and SMM_mean["C2"].iloc[max_idx] < 1) else max_relation)
    data["max_relation_adapted"].append(adapted_max_rel)
    # Add R0 to data: R0 = S(0) * tau_infectious * lambda
    data["R_0"].append(SMM_mean["C1"].iloc[0] *
                       tau_infectious * parameters["lambda"].iloc[0])
    idx_max_mean = SMM_mean["C2"].idxmax()
    data["relation_at_max_mu"].append(np.sqrt(SMM_moments["M020"].iloc[idx_max_mean]) /
                                      SMM_mean["C2"].iloc[idx_max_mean])
    # Add max standard deviation to data
    data["max_std"].append(np.sqrt(SMM_moments["M020"]).max())
    # Add max normalized standard deviation to data
    data["max_std_normalized"].append(np.sqrt(SMM_moments["M020"].max()) / SMM_mean["C2"].max())
    # If ODE model doesn't have last SMM timepoint, skip this time point
    last_timestep = ODE_mean.Time.iloc[-1]
    # Calculate and add mean and variance MSE
    mse_mean = (np.nan if SMM_mean[SMM_mean.Time <= last_timestep]["C2"].isna().any() or ODE_mean["C2"].isna().any()
                else mean_squared_error(SMM_mean[SMM_mean.Time <= last_timestep]["C2"], ODE_mean["C2"]))
    mse_var = (np.nan if SMM_moments[SMM_mean.Time <= last_timestep]["M020"].isna().any() or ODE_moments["M020"].isna().any()
               else mean_squared_error(SMM_moments[SMM_mean.Time <= last_timestep]["M020"], ODE_moments["M020"]))
    data["MSE_mean"].append(mse_mean)
    data["MSE_var"].append(mse_var)
    # Calculate and add mean and variance MAPE
    mape_mean = (np.nan if SMM_mean[SMM_mean.Time <= last_timestep]["C2"].isna().any() or ODE_mean["C2"].isna().any()
                 else mean_absolute_percentage_error(SMM_mean[SMM_mean.Time <= last_timestep]["C2"], ODE_mean["C2"]))
    mape_var = (np.nan if SMM_moments[SMM_mean.Time <= last_timestep]["M020"].isna().any() or ODE_moments["M020"].isna().any()
                else mean_absolute_percentage_error(SMM_moments[SMM_mean.Time <= last_timestep]["M020"], ODE_moments["M020"]))
    data["MAPE_mean"].append(mape_mean)
    data["MAPE_var"].append(mape_var)
    # Calculate and add mean and variance sMAPE
    smape_mean = (np.nan if SMM_mean[SMM_mean.Time <= last_timestep]["C2"].isna().any() or ODE_mean["C2"].isna().any()
                  else scaled_mean_percentage_error(SMM_mean[SMM_mean.Time <= last_timestep]["C2"], ODE_mean["C2"]))
    smape_var = (np.nan if SMM_moments[SMM_mean.Time <= last_timestep]["M020"].isna().any() or ODE_moments["M020"].isna().any()
                 else scaled_mean_percentage_error(SMM_moments[SMM_mean.Time <= last_timestep]["M020"], ODE_moments["M020"]))
    data["sMAPE_mean"].append(smape_mean)
    data["sMAPE_var"].append(smape_var)
    # Calculate and add mean and variance correlation
    corr_mean = SMM_mean[SMM_mean.Time <=
                         last_timestep]["C2"].corr(ODE_mean['C2'])
    corr_var = SMM_moments[SMM_mean.Time <=
                           last_timestep]["M020"].corr(ODE_moments['M020'])
    data["corr_mean"].append(corr_mean)
    data["corr_var"].append(corr_var)

# Write data to file
data_df = pd.DataFrame(data=data)


if not os.path.exists(save_folder):
    os.makedirs(save_folder)

data_df.to_csv(os.path.join(save_folder, "MSEs.csv"), index=False)
