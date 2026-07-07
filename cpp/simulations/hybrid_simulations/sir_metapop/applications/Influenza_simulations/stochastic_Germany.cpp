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

#include "memilio/utils/time_series.h"
#include "simulations/hybrid_simulations/sir_metapop/config/config.cpp"
#include "simulations/hybrid_simulations/sir_metapop/library/smm_helper.h"
#include "memilio/data/analyze_result.h"
#include "memilio/timer/basic_timer.h"
#include "ode_sir/infection_state.h"
#include "smm/simulation.h"
#include "smm/simulation_set.h"
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

template <int n>
void writeMatrixToCSV(const Eigen::Matrix<double, n, n>& mat, const std::string& filename)
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open matrix file\n";
        return;
    }

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            file << mat(i, j);
            if (j < n - 1)
                file << ",";
        }
        file << "\n";
    }

    file.close();
}

template <int n>
std::vector<Eigen::Matrix<double, n, n>>
computeMeanMatrices(std::vector<std::vector<Eigen::Matrix<size_t, n, n>>>& data)
{
    if (data.empty())
        throw std::runtime_error("No data");

    // Check all vectors have same size
    size_t m = data[0].size();
    for (const auto& vec : data) {
        if (vec.size() != m)
            throw std::runtime_error("Vectors have different sizes. Expected: " + std::to_string(m) +
                                     ", got: " + std::to_string(vec.size()));
    }

    size_t numVectors = data.size();

    // Result: one mean matrix per position
    std::vector<Eigen::Matrix<double, n, n>> mean(m, Eigen::Matrix<double, n, n>::Zero());

    // Accumulate
    for (const auto& vec : data) {
        for (size_t k = 0; k < m; ++k) {
            mean[k] += vec[k].template cast<double>();
        }
    }

    // Divide by number of vectors
    for (auto& mat : mean) {
        mat /= static_cast<double>(numVectors);
    }

    return mean;
}

int main()
{
    const size_t num_runs    = 50000;
    const size_t max_order   = 2;
    const auto config        = Config::get_config(Config::ConfigType::ConfigInfluenzaGermany);
    const size_t num_regions = 1;
    if (num_regions != config.num_regions) {
        mio::log_error("Number of regions doesn't match number of regions in config.");
    }

    std::string save_file  = Config::SAVE_DIR + "SMM/";
    auto created_directory = mio::create_directory(save_file);
    save_file += config.name;
    created_directory = mio::create_directory(save_file);
    if (!created_directory) {
        printf("%s\n", created_directory.error().formatted_message().c_str());
        return -1;
    }
    save_file += "/";

    // Initialize model
    auto model = smm_helper::initialize_model<num_regions>(config);
    // Create simulation set
    auto sim_set = mio::smm::SimulationSet<num_regions, mio::osir::InfectionState, max_order>(num_runs, model,
                                                                                              config.t0, config.dt);
    mio::timing::BasicTimer timer;
    timer.start();
    // Advance simulation set
    sim_set.advance(config.tmax);
    timer.stop();

    // Convert result so they fit structure for ensemble_percentile fct
    std::vector<std::vector<mio::TimeSeries<double>>> sim_results;
    auto& all_results = sim_set.get_result();
    //size_t run        = 0;
    for (auto& res : all_results) {
        //(void)res.export_csv(save_file + std::to_string(run) + "_comps.csv");
        sim_results.push_back({res});
        //run += 1;
    }
    // Save percentiles
    auto p05      = mio::ensemble_percentile(sim_results, 0.05);
    auto p25      = mio::ensemble_percentile(sim_results, 0.25);
    auto p50      = mio::ensemble_percentile(sim_results, 0.5);
    auto p75      = mio::ensemble_percentile(sim_results, 0.75);
    auto p95      = mio::ensemble_percentile(sim_results, 0.95);
    auto finished = p05[0].export_csv(save_file + "p05.csv");
    finished      = p25[0].export_csv(save_file + "p25.csv");
    finished      = p50[0].export_csv(save_file + "p50.csv");
    finished      = p75[0].export_csv(save_file + "p75.csv");
    finished      = p95[0].export_csv(save_file + "p95.csv");

    // Save first 100 simulations
    for (size_t sim = 0; sim < 10; ++sim) {
        finished = sim_set.get_result()[sim].export_csv(save_file + std::to_string(sim) + "_result.csv");
    }

    // Save means and moments
    finished = sim_set.get_mean().export_csv(save_file + "means.csv");
    finished = sim_set.get_moments().export_csv(save_file + "moments.csv", sim_set.get_moment_names());
    mio::TimeSeries<double> total_time(1);
    Eigen::VectorXd time = Eigen::VectorXd::Constant(1, mio::timing::time_in_seconds(timer.get_elapsed_time()));
    total_time.add_time_point(0., time);
    auto finished_time = total_time.export_csv(save_file + "total_time.csv", {"Runtime"});

    // Write number new infections to csv file
    std::vector<std::vector<mio::TimeSeries<double>>> new_infections;
    mio::TimeSeries<double> new_infections_mean(num_regions);
    for (auto& sim : sim_set.get_simulations()) {
        mio::TimeSeries<double> new_infections_ts(num_regions);
        for (size_t t = 0; t < sim.num_adoptions.size(); ++t) {
            new_infections_ts.add_time_point(t, Eigen::VectorXd::Zero(num_regions));
            if (static_cast<size_t>(new_infections_mean.get_num_time_points()) <= t) {
                new_infections_mean.add_time_point(t, Eigen::VectorXd::Zero(num_regions));
            }
            for (size_t r = 0; r < num_regions; ++r) {
                double new_inf = sim.num_adoptions[t](r * static_cast<size_t>(mio::osir::InfectionState::Count) +
                                                          static_cast<size_t>(mio::osir::InfectionState::Susceptible),
                                                      r * static_cast<size_t>(mio::osir::InfectionState::Count) +
                                                          static_cast<size_t>(mio::osir::InfectionState::Infected));
                new_infections_ts.get_last_value()[r] = new_inf;
                new_infections_mean.get_value(t)[r] += new_inf;
            }
        }
        new_infections.push_back({new_infections_ts});
    }
    // Normalize mean by number of runs
    for (size_t t = 0; t < static_cast<size_t>(new_infections_mean.get_num_time_points()); ++t) {
        new_infections_mean.get_value(t) /= static_cast<double>(num_runs);
    }
    // Save mean new infections time series
    finished = new_infections_mean.export_csv(save_file + "new_infections_mean.csv");
    // Save new infections percentiles
    auto p05_new_infectuions = mio::ensemble_percentile(new_infections, 0.05);
    auto p25_new_infectuions = mio::ensemble_percentile(new_infections, 0.25);
    auto p50_new_infectuions = mio::ensemble_percentile(new_infections, 0.5);
    auto p75_new_infectuions = mio::ensemble_percentile(new_infections, 0.75);
    auto p95_new_infectuions = mio::ensemble_percentile(new_infections, 0.95);
    finished                 = p05_new_infectuions[0].export_csv(save_file + "p05_new_infections.csv");
    finished                 = p25_new_infectuions[0].export_csv(save_file + "p25_new_infections.csv");
    finished                 = p50_new_infectuions[0].export_csv(save_file + "p50_new_infections.csv");
    finished                 = p75_new_infectuions[0].export_csv(save_file + "p75_new_infections.csv");
    finished                 = p95_new_infectuions[0].export_csv(save_file + "p95_new_infections.csv");

    std::cout << "SMM Elapsed time: " << timer.get_elapsed_time() << std::endl << std::flush;

    return 0;
}
