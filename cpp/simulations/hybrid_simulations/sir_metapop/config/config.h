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
const std::string SAVE_DIR = "/hpc_data/bick_ju/TemporalHybrid/";

const std::vector<std::string> closure_string = {"truncation", "pair_approx", "lognorm"};

struct Config {
    // Config name
    std::string name;
    size_t num_regions;
    double t0, tmax, dt;
    // Lambdas (S->E) for every region
    std::vector<double> lambdas;
    // Recovery rate (I->R) is the same for all regions
    double gamma;
    // Number of initially infected for given regions. If a region id is not present in the vector, the number of initially infected in that region is 0.
    std::vector<std::pair<int, double>> I0s;
    // Total population for every region
    std::vector<double> total_populations;
    // Transition rates
    std::vector<mio::smm::TransitionRate<ScalarType, mio::osir::InfectionState>> transition_rates;
};

enum class ConfigType
{
    Config1,
    Config2,
    Config3,
    Config4,
    Config3r1,
    Config3r2,
    Config3r3
};

Config get_config(ConfigType type);

} // namespace Config

#endif // CONFIG_H
