/* 
* Copyright (C) 2020-2026 MEmilio
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

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include <cstddef>
#include <vector>

// Macro to skip main() in stochastic_Germany.cpp
#define STOCHASTIC_SIM_BINDINGS_SKIP_MAIN
#include "simulations/hybrid_simulations/sir_metapop/applications/Influenza_simulations/stochastic_Germany.cpp"
#include "simulations/hybrid_simulations/sir_metapop/config/config.h"

namespace py = pybind11;

std::vector<std::vector<double>> run_germany(int full_or_scaled, bool use_new_inf, size_t num_runs, double lambda,
                                             double gamma, double nu, double I0, double R0, std::vector<double> peaks,
                                             std::vector<double> rhos, std::vector<double> sigmas)
{
    Config::Config config;
    if (full_or_scaled) {
        config = Config::get_config(Config::ConfigType::ConfigInfluenzaGermany_full);
    }
    else {
        config = Config::get_config(Config::ConfigType::ConfigInfluenzaGermany);
    }
    config.lambdas[0]         = lambda;
    config.gamma              = gamma;
    config.nu                 = nu;
    config.I0s[0]             = {0, I0};
    config.R0s[0]             = {0, R0};
    config.season_peaks       = peaks;
    config.seasonality_rhos   = rhos;
    config.seasonality_sigmas = sigmas;
    return run_stochastic_simulation_set_one_region(num_runs, config, false, use_new_inf);
}

std::vector<std::vector<double>> run_age_groups(int full_or_scaled, bool use_new_inf, size_t num_runs,
                                                std::vector<double> lambdas, double gamma, double nu,
                                                std::vector<double> I0s, std::vector<double> R0s,
                                                std::vector<double> peaks, std::vector<double> rhos,
                                                std::vector<double> sigmas)
{
    Config::Config config;
    if (full_or_scaled) {
        config = Config::get_config(Config::ConfigType::ConfigInfluenzaAgeGroups_full);
    }
    else {
        config = Config::get_config(Config::ConfigType::ConfigInfluenzaAgeGroups);
    }
    config.lambdas = lambdas;
    config.gamma   = gamma;
    config.nu      = nu;
    for (size_t ag = 0; ag < I0s.size(); ++ag) {
        config.I0s[ag] = {ag, I0s[ag]};
        config.R0s[ag] = {ag, R0s[ag]};
    }
    config.season_peaks       = peaks;
    config.seasonality_rhos   = rhos;
    config.seasonality_sigmas = sigmas;
    return run_stochastic_simulation_set_five_regions(num_runs, config, false, use_new_inf);
}

std::vector<std::vector<double>> run_regions(int full_or_scaled, bool use_new_inf, size_t num_runs,
                                             std::vector<double> lambdas, double gamma, double nu,
                                             std::vector<double> I0s, std::vector<double> R0s,
                                             std::vector<double> peaks, std::vector<double> rhos,
                                             std::vector<double> sigmas)
{
    Config::Config config;
    if (full_or_scaled) {
        config = Config::get_config(Config::ConfigType::ConfigInfluenzaRegions_full);
    }
    else {
        config = Config::get_config(Config::ConfigType::ConfigInfluenzaRegions);
    }
    config.lambdas = lambdas;
    config.gamma   = gamma;
    config.nu      = nu;
    for (size_t r = 0; r < I0s.size(); ++r) {
        config.I0s[r] = {r, I0s[r]};
        config.R0s[r] = {r, R0s[r]};
    }
    config.season_peaks       = peaks;
    config.seasonality_rhos   = rhos;
    config.seasonality_sigmas = sigmas;
    return run_stochastic_simulation_set_four_regions(num_runs, config, false, use_new_inf);
}

PYBIND11_MODULE(_simulation_stochastic, m)
{
    m.def("run_germany", &run_germany, "Simulate simulation set for Germany", py::arg("full_or_scaled"),
          py::arg("use_new_inf"), py::arg("num_runs"), py::arg("lambda"), py::arg("gamma"), py::arg("nu"),
          py::arg("I0"), py::arg("R0"), py::arg("peaks"), py::arg("rhos"), py::arg("sigmas"));
    m.def("run_age_groups", &run_age_groups, "Simulate simulation set for Germany age-resolved",
          py::arg("full_or_scaled"), py::arg("use_new_inf"), py::arg("num_runs"), py::arg("lambdas"), py::arg("gamma"),
          py::arg("nu"), py::arg("I0s"), py::arg("R0s"), py::arg("peaks"), py::arg("rhos"), py::arg("sigmas"));
    m.def("run_regions", &run_regions, "Simulate simulation set for Germany region-resolved", py::arg("full_or_scaled"),
          py::arg("use_new_inf"), py::arg("num_runs"), py::arg("lambdas"), py::arg("gamma"), py::arg("nu"),
          py::arg("I0s"), py::arg("R0s"), py::arg("peaks"), py::arg("rhos"), py::arg("sigmas"));
    m.attr("__version__") = "dev";
}
