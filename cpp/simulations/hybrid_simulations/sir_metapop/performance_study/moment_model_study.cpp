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

#include "memilio/timer/basic_timer.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "simulations/hybrid_simulations/sir_metapop/config/config.cpp"
#include "simulations/hybrid_simulations/sir_metapop/config/config.h"
#include "smm_moments/model.h"
#include "smm_moments/simulation.h"
#include "simulations/hybrid_simulations/sir_metapop/library/moment_helper.h"
#include "memilio/data/analyze_result.h"
#include "smm_moments/closure_functions.h"
#include <cstddef>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

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

std::vector<double> read_parameter_csv(std::string filename)
{
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Failed to open parameter file\n";
    }

    std::string line;

    // Skip header
    std::getline(file, line);

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string I0_str, R0_str, lambda_str;

        std::getline(ss, I0_str, ',');
        std::getline(ss, R0_str, ',');
        std::getline(ss, lambda_str, ',');

        double I0     = std::stod(I0_str);
        double R0     = std::stod(R0_str);
        double lambda = std::stod(lambda_str);
        return {I0, R0, lambda};
    }
    return {};
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
    auto config                = Config::get_config(Config::ConfigType::ConfigPerformanceStudySIR);
    const size_t num_regions   = 1;
    const size_t closure_order = 3;
    double min_step_size       = 0.0001;
    double init_time           = 0.0;
    auto closure_func          = &mio::smm_moments::truncation_closure<num_regions, closure_order>;
    if (num_regions != config.num_regions) {
        mio::log_error("Number of regions doesn't match number of regions in config.");
    }

    std::vector<size_t> I0_boundaries     = {0, static_cast<size_t>(0.01 * config.total_populations[0])};
    std::vector<size_t> R0_boundaries     = {0, static_cast<size_t>(0.7 * config.total_populations[0])};
    std::vector<double> lambda_boundaries = {0.0000014, 0.000006};

    size_t num_samples     = 1000;
    std::string save_file  = Config::SAVE_DIR + "Moments/";
    auto created_directory = mio::create_directory(save_file);
    save_file += config.name;
    created_directory = mio::create_directory(save_file);
    if (!created_directory) {
        printf("%s\n", created_directory.error().formatted_message().c_str());
        return -1;
    }
    save_file += "/";

    std::string init_dir_base = "/Users/julia/sim_outputs/output/SMM/performance_study_SIR/";

    for (size_t sample = 0; sample < num_samples; ++sample) {
        // Create sample directory
        std::string save_file_sample = save_file + "sample_" + std::to_string(sample);
        created_directory            = mio::create_directory(save_file_sample);
        if (!created_directory) {
            printf("%s\n", created_directory.error().formatted_message().c_str());
            return -1;
        }
        save_file_sample += "/";

        std::string init_dir = init_dir_base + "sample_" + std::to_string(sample) + "/";

        std::string file_expected_values = init_dir + "means.csv";
        std::string file_moment_values   = init_dir + "moments.csv";

        auto params = read_parameter_csv(init_dir + "parameters.csv");

        // Sample I0, R0 and lambda
        config.I0s[0].second = params[0];
        config.R0s[0].second = params[1];
        config.lambdas[0]    = params[2];

        write_parameter_csv(config.I0s[0].second, config.R0s[0].second, config.lambdas[0],
                            save_file_sample + "parameters.csv");

        // Read init expected values and moments
        auto expected_values_init = moment_helper::read_expected_values(file_expected_values, init_time);
        auto moments_init = moment_helper::read_moments<num_regions, closure_order>(file_moment_values, init_time);

        // Initialize model
        auto model = moment_helper::initialize_model<num_regions, closure_order>(expected_values_init, moments_init,
                                                                                 config, closure_func);

        // Create simulation
        auto sim = mio::smm_moments::Simulation<num_regions, closure_order>(model, init_time, config.dt);

        sim.get_integrator_core().get_dt_max() = config.dt;
        if (min_step_size > 0) {
            sim.get_integrator_core().get_dt_min() = min_step_size;
        }

        mio::timing::BasicTimer timer;
        timer.start();
        // Advance simulation set
        sim.advance(config.tmax);
        timer.stop();

        // Save means and moments
        auto means   = sim.get_expected_values_time_series();
        auto moments = sim.get_moment_time_series(closure_order);

        int num_steps = static_cast<int>((config.tmax - init_time) / config.dt) + 1;
        std::vector<double> interpolation_tps(num_steps);
        for (int i = 0; i < num_steps; ++i) {
            interpolation_tps[i] = init_time + i * config.dt;
        }
        means          = mio::interpolate_simulation_result(means, interpolation_tps);
        auto done      = means.export_csv(save_file_sample + "means.csv");
        auto moment_ts = mio::interpolate_simulation_result(moments.first, interpolation_tps);
        done           = moment_ts.export_csv(save_file_sample + "moments.csv", moments.second);
        // Save time
        write_time_csv(mio::timing::time_in_seconds(timer.get_elapsed_time()), save_file_sample + "total_time.csv");
    }

    return 0;
}
