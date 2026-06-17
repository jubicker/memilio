/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Julia Bicker
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/
#ifndef CONFIG_H
#define CONFIG_H

#include "memilio/config.h"
#include "ode_sir/infection_state.h"
#include "smm/parameters.h"
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace Config
{
const std::string SAVE_DIR =
    "/p/project1/loki/bicker1/memilio/output/"; //"/Users/julia/repos/fork/memilio/output/";//"/hpc_data/bick_ju/TemporalHybrid/";

const std::vector<std::string> closure_string          = {"truncation", "pair_approx", "lognorm", "lognorm_zero_infl"};
const std::vector<std::string> switch_condition_string = {"mean_threshold_rel", "var_gradient", "mean_threshold_abs",
                                                          "mean_stddev_relation"};

struct Config {
    // Config name
    std::string name;
    size_t num_regions;
    double t0, tmax, dt;
    // Lambdas (S->E) for every region
    std::vector<double> lambdas;
    // Recovery rate (I->R) is the same for all regions
    double gamma;
    // Immunity loss rate (R->S) is the same for all regions
    double nu;
    // Number of initially infected for given regions. If a region id is not present in the vector, the number of initially infected in that region is 0.
    std::vector<std::pair<int, double>> I0s;
    // Number of initially recovered for given regions. If a region id is not present in the vector, the number of initially recovered in that region is 0.
    std::vector<std::pair<int, double>> R0s;
    // Total population for every region
    std::vector<double> total_populations;
    // Transition rates
    std::vector<mio::smm::TransitionRate<ScalarType, mio::osir::InfectionState>> transition_rates;
    // Seasonality parameters
    // Start day of the first season given in days from 1st January in the first season, e.g. 0 for 1st January, -185 for 1st July (in the year before 1st January of that season), etc.
    int first_season_start_day;
    // Time points of peak transmission rate for each season given in days from 1st January in the first season, e.g. 0 for 1st January, 10 for 11th January, etc.
    std::vector<int> season_peaks;
    // Variation factor for each season
    std::vector<double> seasonality_rhos;
    // Standard deviation for the Gaussian function for each season
    std::vector<double> seasonality_sigmas;
};

enum class ConfigType
{
    ConfigSIRVaryI0NoExchange,
    ConfigSIRVaryI0Exchange,
    ConfigSIRVaryR0NoExchange,
    ConfigSIRVaryR0Exchange,
    ConfigSIRSVaryI0NoExchange,
    ConfigSIRSVaryI0Exchange,
    ConfigSIRSVaryR0NoExchange,
    ConfigSIRSVaryR0Exchange,
    ConfigSIRSVaryI0NoExchange_Seasonal,
    ConfigSIRSVaryI0Exchange_Seasonal,
    ConfigSIRSVaryR0NoExchange_Seasonal,
    ConfigSIRSVaryR0Exchange_Seasonal,
    ConfigTestExchange,
    ConfigTest2ndOutbreak,
    ConfigPerformanceStudySIR,
    ConfigPerformanceStudySIRS,
    ConfigPerformanceStudySIRS_seasonal,
};

Config get_config(ConfigType type);

} // namespace Config

#endif // CONFIG_H
