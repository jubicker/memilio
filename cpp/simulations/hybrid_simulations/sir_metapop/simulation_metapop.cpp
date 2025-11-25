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

#include "config/config.cpp"
#include "library/smm_helper.h"
#include "memilio/data/analyze_result.h"
#include "memilio/timer/basic_timer.h"
#include "ode_sir/infection_state.h"
#include "smm/simulation.h"
#include <cstddef>
#include <string>
#include <vector>
#include <omp.h>

template <int NumRegions>
mio::TimeSeries<double> run_smm_sim(int sim_num, std::string save_file, const Config::Config& config)
{
    auto model = smm_helper::initialize_model<NumRegions>(config);
    auto sim   = mio::smm::Simulation(model, config.t0, config.dt);
    sim.advance(config.tmax);

    std::string output_file = save_file + std::to_string(sim_num) + "_comps.csv";

    int num_steps = static_cast<int>(config.tmax / config.dt) + 1;

    std::vector<double> interpolation_tps(num_steps);

    for (int i = 0; i < num_steps; ++i) {
        interpolation_tps[i] = i * config.dt;
    }
    auto result = mio::interpolate_simulation_result(sim.get_result(), interpolation_tps);
    std::vector<std::string> names(NumRegions * 3);
    for (int r = 0; r < NumRegions; ++r) {
        names.push_back("S" + std::to_string(r));
        names.push_back("I" + std::to_string(r));
        names.push_back("R" + std::to_string(r));
    }
    auto done = result.export_csv(output_file, names);
    return result;
}

int main()
{
    const size_t num_runs    = 1000;
    const size_t max_order   = 3;
    const auto config        = Config::get_config(Config::ConfigType::Config1);
    const size_t num_regions = 1;
    std::vector<mio::TimeSeries<double>> sim_results(
        num_runs, mio::TimeSeries<double>(static_cast<size_t>(mio::osir::InfectionState::Count) * num_regions));
    std::string save_file = Config::SAVE_DIR + "SMM/";
    save_file += config.name;
    auto created_directory = mio::create_directory(save_file);

    if (!created_directory) {
        printf("%s\n", created_directory.error().formatted_message().c_str());
        return -1;
    }

    std::vector<double> time(num_runs);
    save_file += "/";

// Run multiple simulations
#pragma omp parallel for
    for (size_t run = 0; run < num_runs; ++run) {
        mio::timing::BasicTimer timer;
        timer.start();
        sim_results[run] = run_smm_sim<num_regions>(run, save_file, config);
        timer.stop();
        time[run] = timer.get_elapsed_time();
    }
#pragma omp single
    {
        // Convert time to timeseries to make exporting to csv easier
        mio::TimeSeries<double> time_ts(1);
        for (size_t i = 0; i < num_runs; i++) {
            time_ts.add_time_point(i, Eigen::VectorXd::Constant(1, time[i]));
        }
        auto finished_time = time_ts.export_csv(save_file + "runtimes.csv", {"runtime"});
    }
    // Calculate moments and expected values from the simulation results
    auto means       = smm_helper::calculate_means_from_sim<num_regions>(sim_results);
    auto mean_string = means.second;
    auto moments     = smm_helper::calculate_moments_from_sim<max_order, num_regions>(sim_results);
    auto finished    = means.first.export_csv(save_file + "means.csv", means.second);
    finished         = moments.first.export_csv(save_file + "moments.csv", moments.second);
}
