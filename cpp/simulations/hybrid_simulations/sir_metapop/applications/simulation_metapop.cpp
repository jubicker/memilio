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

#include "simulations/hybrid_simulations/sir_metapop/config/config.cpp"
#include "simulations/hybrid_simulations/sir_metapop/library/smm_helper.h"
#include "memilio/data/analyze_result.h"
#include "memilio/timer/basic_timer.h"
#include "ode_sir/infection_state.h"
#include "smm/simulation.h"
#include "smm/simulation_set.h"
#include <cstddef>
#include <string>
#include <vector>
#include <omp.h>

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
    const auto config        = Config::get_config(Config::ConfigType::ConfigSIRVaryI0Exchange);
    const size_t num_regions = 4;
    if (num_regions != config.num_regions) {
        mio::log_error("Number of regions doesn't match number of regions in config.");
    }

    std::string save_file = Config::SAVE_DIR + "SMM/";
    save_file += config.name;
    auto created_directory = mio::create_directory(save_file);
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
    auto p00      = mio::ensemble_percentile(sim_results, 0.0);
    auto p05      = mio::ensemble_percentile(sim_results, 0.05);
    auto p25      = mio::ensemble_percentile(sim_results, 0.25);
    auto p50      = mio::ensemble_percentile(sim_results, 0.5);
    auto p75      = mio::ensemble_percentile(sim_results, 0.75);
    auto p95      = mio::ensemble_percentile(sim_results, 0.95);
    auto p100     = mio::ensemble_percentile(sim_results, 1.0);
    auto finished = p00[0].export_csv(save_file + "p00.csv");
    finished      = p05[0].export_csv(save_file + "p05.csv");
    finished      = p25[0].export_csv(save_file + "p25.csv");
    finished      = p50[0].export_csv(save_file + "p50.csv");
    finished      = p75[0].export_csv(save_file + "p75.csv");
    finished      = p95[0].export_csv(save_file + "p95.csv");
    finished      = p100[0].export_csv(save_file + "p100.csv");

    // Save first 100 simulations
    for (size_t sim = 0; sim < 100; ++sim) {
        finished = sim_set.get_result()[sim].export_csv(save_file + std::to_string(sim) + "_result.csv");
    }

    // Save means and moments
    finished = sim_set.get_mean().export_csv(save_file + "means.csv");
    finished = sim_set.get_moments().export_csv(save_file + "moments.csv", sim_set.get_moment_names());
    // Save times
    mio::TimeSeries<double> time_ts(1);
    auto& all_times = sim_set.get_sim_times();
    for (size_t i = 0; i < all_times.size(); i++) {
        Eigen::VectorXd time = Eigen::VectorXd::Zero(1);
        time[0]              = all_times[i];
        time_ts.add_time_point(i, time);
    }
    auto finished_time = time_ts.export_csv(save_file + "runtimes.csv", {"Runtime"});
    mio::TimeSeries<double> total_time(1);
    Eigen::VectorXd time = Eigen::VectorXd::Constant(1, mio::timing::time_in_seconds(timer.get_elapsed_time()));
    total_time.add_time_point(0., time);
    finished_time = total_time.export_csv(save_file + "total_time.csv", {"Runtime"});

    // TODO - Write number of transitions to csv file
    if (false) {
        std::vector<
            std::vector<Eigen::Matrix<size_t, static_cast<size_t>(mio::osir::InfectionState::Count) * num_regions,
                                      static_cast<size_t>(mio::osir::InfectionState::Count) * num_regions>>>
            number_transitions;
        number_transitions.reserve(num_runs);
        for (auto& sim : sim_set.get_simulations()) {
            number_transitions.push_back(sim.num_transitions);
        }
        auto mean = computeMeanMatrices(number_transitions);
        for (size_t t = 0; t < mean.size(); ++t) {
            writeMatrixToCSV(mean[t], save_file + std::to_string(t) + "_transitions.csv");
        }
    }

    // TODO - Write number of total events to csv file
    // Convert result so they fit structure for ensemble_percentile fct
    std::vector<std::vector<mio::TimeSeries<double>>> events;
    for (size_t sim = 0; sim < sim_set.get_simulations().size(); ++sim) {
        events.push_back({sim_set.get_simulations()[sim].num_events});
    }
    // Save percentiles
    auto p50_events  = mio::ensemble_percentile(events, 0.5);
    auto p00_events  = mio::ensemble_percentile(events, 0.0);
    auto p100_events = mio::ensemble_percentile(events, 1.0);
    finished         = p50_events[0].export_csv(save_file + "events_median.csv");
    finished         = p00_events[0].export_csv(save_file + "events_min.csv");
    finished         = p100_events[0].export_csv(save_file + "events_max.csv");

    std::cout << "SMM Elapsed time: " << timer.get_elapsed_time() << std::endl << std::flush;

    return 0;
}
