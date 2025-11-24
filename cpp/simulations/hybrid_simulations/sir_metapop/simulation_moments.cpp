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

#include "hybrid_simulations/sir_metapop/config/config.cpp"
#include "hybrid_simulations/sir_metapop/library/moments/model.h"
#include "hybrid_simulations/sir_metapop/library/moments/simulation.h"
#include "hybrid_simulations/sir_metapop/library/moment_helper.h"
#include <cstddef>
#include <string>

template <size_t NumRegions, size_t ClosureOrder>
void run_moments_simulation(std::string save_dir, const Config::Config& config,
                            Eigen::Array<double, Eigen::Dynamic, 1>& expected_values_init,
                            Eigen::Array<double, Eigen::Dynamic, 1>& moments_init, double min_step_size)
{
    auto model = moment_helper::initialize_model<NumRegions, ClosureOrder>(expected_values_init, moments_init, config);
    auto sim   = smm_moments::Simulation<NumRegions, ClosureOrder>(model, config.t0, config.dt);
    if (min_step_size > 0) {
        sim.get_integrator().get_dt_min() = min_step_size;
    }
    sim.advance(config.tmax);

    //Save expected values time series
    std::string output_file_expected_values = save_dir + "/expected_values.csv";
    auto done = sim.get_expected_values_time_series().export_csv(output_file_expected_values);
    // Save moment time series
    auto moments                    = sim.get_moment_time_series(ClosureOrder);
    std::string output_file_moments = save_dir + "/moments.csv";
    done                            = moments.first.export_csv(output_file_moments, moments.second);
}

int main()
{
    auto config                = Config::get_config(Config::ConfigType::Config1);
    const size_t closure_order = 3;
    const size_t num_regions   = 1;
    double min_step_size       = 0.00001;
    double init_time           = 0.0;
    std::string save_file      = Config::SAVE_DIR + "Moments/";
    save_file += config.name; //+ "/closure_order_" + std::to_string(closure_order) + "/" + std::to_string(init_time);
    auto created_directory = mio::create_directory(save_file);
    if (!created_directory) {
        printf("%s\n", created_directory.error().formatted_message().c_str());
        return -1;
    }
    save_file += "/closure_order_" + std::to_string(closure_order);
    created_directory = mio::create_directory(save_file);
    if (!created_directory) {
        printf("%s\n", created_directory.error().formatted_message().c_str());
        return -1;
    }
    save_file += "/" + std::to_string(init_time);
    created_directory = mio::create_directory(save_file);
    if (!created_directory) {
        printf("%s\n", created_directory.error().formatted_message().c_str());
        return -1;
    }

    config.t0 -= init_time;
    std::string file_expected_values = Config::SAVE_DIR + "SMM/" + config.name + "/means.csv";
    std::string file_moment_values   = Config::SAVE_DIR + "SMM/" + config.name + "/moments.csv";
    auto ecpected_values             = moment_helper::read_expected_values(file_expected_values, init_time);
    auto moments = moment_helper::read_moments<num_regions, closure_order>(file_moment_values, init_time);
    run_moments_simulation<num_regions, closure_order>(save_file, config, ecpected_values, moments, min_step_size);
    return 0;
}
