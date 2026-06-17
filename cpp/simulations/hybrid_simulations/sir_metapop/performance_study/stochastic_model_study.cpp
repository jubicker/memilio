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
#include <fstream>
#include <vector>
#include <omp.h>

void write_parameter_csv(double I0, double R0, double lambda, std::string filename)
{
    std::ofstream file(filename);

    if (!file.is_open())
        return;

    // Header
    file << "I0,R0,lambda\n";

    // Data row
    file << I0 << "," << R0 << "," << lambda << "\n";

    file.close();
}

void write_time_csv(double time, std::string filename)
{
    std::ofstream file(filename);

    if (!file.is_open())
        return;

    // Header
    file << "Time,Runtime\n";

    // Data row
    file << "0.0," << time << "\n";

    file.close();
}

int main()
{
    const size_t max_order   = 2;
    const size_t num_runs    = 50000;
    auto config              = Config::get_config(Config::ConfigType::ConfigPerformanceStudySIR);
    const size_t num_regions = 1;
    if (num_regions != config.num_regions) {
        mio::log_error("Number of regions doesn't match number of regions in config.");
    }

    std::vector<size_t> I0_boundaries     = {0, static_cast<size_t>(0.01 * config.total_populations[0])};
    std::vector<size_t> R0_boundaries     = {0, static_cast<size_t>(0.7 * config.total_populations[0])};
    std::vector<double> lambda_boundaries = {0.0000014, 0.000006};

    size_t num_samples     = 1000;
    std::string save_file  = Config::SAVE_DIR + "SMM/";
    auto created_directory = mio::create_directory(save_file);
    save_file += config.name;
    created_directory = mio::create_directory(save_file);
    if (!created_directory) {
        printf("%s\n", created_directory.error().formatted_message().c_str());
        return -1;
    }
    save_file += "/";

    for (size_t sample = 0; sample < num_samples; ++sample) {
        // Create sample directory
        std::string save_file_sample = save_file + "sample_" + std::to_string(sample);
        created_directory            = mio::create_directory(save_file_sample);
        if (!created_directory) {
            printf("%s\n", created_directory.error().formatted_message().c_str());
            return -1;
        }
        save_file_sample += "/";

        // Sample I0, R0 and lambda
        config.I0s[0].second = mio::UniformIntDistribution<size_t>::get_instance()(mio::thread_local_rng(),
                                                                                   I0_boundaries[0], I0_boundaries[1]);
        config.R0s[0].second = mio::UniformIntDistribution<size_t>::get_instance()(mio::thread_local_rng(),
                                                                                   R0_boundaries[0], R0_boundaries[1]);
        config.lambdas[0]    = mio::UniformDistribution<ScalarType>::get_instance()(
            mio::thread_local_rng(), lambda_boundaries[0], lambda_boundaries[1]);

        write_parameter_csv(config.I0s[0].second, config.R0s[0].second, config.lambdas[0],
                            save_file_sample + "parameters.csv");

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

        // Save means and moments
        auto finished = sim_set.get_mean().export_csv(save_file_sample + "means.csv");
        finished      = sim_set.get_moments().export_csv(save_file_sample + "moments.csv", sim_set.get_moment_names());
        // Save time
        write_time_csv(mio::timing::time_in_seconds(timer.get_elapsed_time()), save_file_sample + "total_time.csv");
    }

    return 0;
}
