import pandas as pd
import numpy as np
from sklearn.metrics import mean_squared_error
from sklearn.metrics import mean_absolute_percentage_error
import os


def scaled_mean_percentage_error(y_true, y_pred):
    denominator = (np.abs(y_true) + np.abs(y_pred)) / 2.
    nonzero = denominator != 0
    if not nonzero.any():
        return np.nan
    return np.mean(np.abs(y_pred[nonzero] - y_true[nonzero]) / denominator[nonzero])


sample_start = 0
sample_end = 4319

data = {"I_init_mean": [], "R_init_mean": [], "lambda": [], "R0_mean": [], "std_I_init": [],
        "std_R_init": [], "std_I_init_norm": [], "std_R_init_norm": [], "sMAPE_mean": [], "sMAPE_var": [],
        "MSE_mean": [], "MSE_var": [], "NaN_mean": [], "NaN_var": []}

base_folder = "/Users/julia/sim_outputs/output/"
smm_folder = base_folder + "SMM/performance_study_SIR/"
moment_folder = base_folder + "Moments/performance_study_SIR"
save_folder = base_folder + "ParameterStudy/SIR/"

subfolder = [0, 10, 20, 40, 60, 80, 100]

for f in subfolder:
    for s in range(sample_start, sample_end + 1):
        moment_subfolder = moment_folder + f"{f}/"
        # Read in SMM mean and moments
        SMM_mean = pd.read_csv(smm_folder + "sample_" + str(s) + "/means.csv")
        SMM_moments = pd.read_csv(
            smm_folder + "sample_" + str(s) + "/moments.csv")
        # Read in ODE mean and moments
        ODE_mean = pd.read_csv(
            moment_subfolder + "sample_" + str(s) + "/means.csv")
        ODE_moments = pd.read_csv(
            moment_subfolder + "sample_" + str(s) + "/moments.csv")
        # Read in parameters
        parameters = pd.read_csv(moment_subfolder + "sample_" +
                                 str(s) + "/parameters.csv")
        # Add parameters to data
        data["I_init_mean"].append(parameters["I_init_mean"].iloc[0])
        data["R_init_mean"].append(parameters["R_init_mean"].iloc[0])
        data["lambda"].append(parameters["lambda"].iloc[0])
        data["R0_mean"].append(parameters["R0_mean"].iloc[0])
        data["std_I_init"].append(parameters["std_I_init"].iloc[0])
        data["std_R_init"].append(parameters[" std_R_init"].iloc[0])
        # Add normalized parameters
        data["std_I_init_norm"].append(
            parameters["std_I_init"].iloc[0]/parameters["I_init_mean"].iloc[0])
        data["std_R_init_norm"].append(
            parameters[" std_R_init"].iloc[0]/parameters["R_init_mean"].iloc[0])
        # Get first an last time point of ODE model
        first_tp = ODE_mean.Time.iloc[0]
        last_tp = ODE_mean.Time.iloc[-1]
        # Add NaN value or sMAPE/MSE for mean
        if (ODE_mean.C2.isna().any()):
            data["sMAPE_mean"].append(np.nan)
            data["MSE_mean"].append(np.nan)
            data["NaN_mean"].append(True)
        else:
            smm_values = SMM_mean[(SMM_mean.Time >= first_tp) & (
                SMM_mean.Time <= last_tp)].C2
            ode_values = ODE_mean[(ODE_mean.Time >= first_tp) & (
                ODE_mean.Time <= last_tp)].C2
            error = scaled_mean_percentage_error(
                np.array(smm_values), np.array(ode_values))
            data["sMAPE_mean"].append(error)
            mse = mean_squared_error(
                np.array(ode_values), np.array(smm_values))
            data["MSE_mean"].append(mse)
            data["NaN_mean"].append(False)
        # Add NaN value or sMAPE/MSE for var
        if (ODE_moments["M020"].isna().any()):
            data["sMAPE_var"].append(np.nan)
            data["MSE_var"].append(np.nan)
            data["NaN_var"].append(True)
        else:
            smm_values = SMM_moments[(SMM_moments.Time >= first_tp) & (
                SMM_moments.Time <= last_tp)]["M020"]
            ode_values = ODE_moments[(ODE_moments.Time >= first_tp) & (
                ODE_moments.Time <= last_tp)]["M020"]
            error = scaled_mean_percentage_error(
                np.array(smm_values), np.array(ode_values))
            data["sMAPE_var"].append(error)
            mse = mean_squared_error(
                np.array(ode_values), np.array(smm_values))
            data["MSE_var"].append(mse)
            data["NaN_var"].append(False)


# Write data to file
data_df = pd.DataFrame(data=data)


if not os.path.exists(save_folder):
    os.makedirs(save_folder)

data_df.to_csv(os.path.join(save_folder, "results.csv"), index=False)
