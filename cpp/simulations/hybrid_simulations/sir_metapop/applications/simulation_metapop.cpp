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

/**
 * @brief Run one simulation of the smm. The comps.csv is saved as result.
 * @tparam NumRegions Number of regions.
 * @param[in] sim_num Simulation number used as seed for the simulation.
 * @param[in] save_file File the output time series is saved to.
 * @param[in] config Config used for the simulation i.e. parameters and initial populations.
 * @param[in, out] timer_init Timer used to time initialization.
 * @param[in, out] timer_sim Timer used to time simulation.
 */
template <int NumRegions>
mio::TimeSeries<double> run_smm_sim(int sim_num, std::string save_file, const Config::Config& config,
                                    mio::timing::BasicTimer& timer_init, mio::timing::BasicTimer& timer_sim)
{
    timer_init.start();
    // Initialize model
    auto model = smm_helper::initialize_model<NumRegions>(config);
    // Create simulation
    auto sim = mio::smm::Simulation(model, config.t0, config.dt);
    timer_init.stop();
    timer_sim.start();
    // Advance simulation until tmax
    sim.advance(config.tmax);
    timer_sim.stop();

    std::string output_file = save_file + std::to_string(sim_num) + "_comps.csv";

    // Output is interpolated to time steps of size dt
    int num_steps = static_cast<int>(config.tmax / config.dt) + 1;
    std::vector<double> interpolation_tps(num_steps);

    for (int i = 0; i < num_steps; ++i) {
        interpolation_tps[i] = i * config.dt;
    }
    auto result = mio::interpolate_simulation_result(sim.get_result(), interpolation_tps);
    auto done   = result.export_csv(output_file);
    return result;
}

int main()
{
    const size_t num_runs    = 10000;
    const size_t max_order   = 3;
    const auto config        = Config::get_config(Config::ConfigType::Config1);
    const size_t num_regions = 1;
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
    // Advance simulation set
    sim_set.advance(config.tmax);
    sim_set.calculate_outputs();

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
    Eigen::VectorXd time = Eigen::VectorXd::Constant(1, sim_set.get_total_time());
    total_time.add_time_point(0., time);
    finished_time = total_time.export_csv(save_file + "total_time.csv", {"Runtime"});

    return 0;
}
