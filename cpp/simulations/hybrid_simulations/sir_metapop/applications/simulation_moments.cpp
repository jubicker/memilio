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

int main()
{
    auto config                = Config::sir::get_config(Config::sir::ConfigType::ConfigConference1);
    const size_t closure_order = 2;
    const size_t num_regions   = 1;
    double min_step_size       = 0.0001;
    double init_time           = 0.0;
    size_t closure             = 0;
    auto closure_func          = &mio::smm_moments::truncation_closure<num_regions, closure_order>;
    if (closure == 0) {
        closure_func = &mio::smm_moments::truncation_closure<num_regions, closure_order>;
    }
    else if (closure == 1) {
        closure_func = &mio::smm_moments::pairapprox_closure<num_regions, closure_order>;
    }
    else if (closure == 2) {
        closure_func = &mio::smm_moments::lognormal_closure<num_regions, closure_order>;
    }
    else if (closure == 3) {
        closure_func = &mio::smm_moments::lognormal_zero_inflation_closure<num_regions, closure_order>;
    }
    else {
        mio::log_error("Unkown closure type: ", closure);
    }

    if (num_regions != config.num_regions) {
        mio::log_error("Number of regions doesn't match number of regions in config.");
    }

    // INPUT - SIR/SIRS
    std::string save_file = Config::SAVE_DIR + "conference/ODE/";
    save_file += config.name;
    auto created_directory = mio::create_directory(save_file);
    if (!created_directory) {
        printf("%s\n", created_directory.error().formatted_message().c_str());
        return -1;
    }
    save_file += "/" + Config::closure_string[closure];
    created_directory = mio::create_directory(save_file);
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
    save_file += "/";

    // Start at init time
    config.t0 += init_time;

    // INPUT - SIR/SIRS
    std::string init_dir = Config::SAVE_DIR + "conference/SMM/" + config.name + "/";

    std::string file_expected_values = init_dir + "means.csv";
    std::string file_moment_values   = init_dir + "moments.csv";

    // Read init expected values and moments
    auto expected_values_init = moment_helper::read_expected_values(file_expected_values, init_time);
    auto moments_init         = moment_helper::read_moments<num_regions, closure_order>(file_moment_values, init_time);

    // Initialize model
    auto model = moment_helper::initialize_model<num_regions, closure_order>(expected_values_init, moments_init, config,
                                                                             closure_func);
    // Create simulation
    auto sim = mio::smm_moments::Simulation<num_regions, closure_order>(model, init_time, config.dt);

    sim.get_integrator_core().get_dt_max() = config.dt;
    if (min_step_size > 0) {
        sim.get_integrator_core().get_dt_min() = min_step_size;
    }

    mio::timing::BasicTimer timer;
    timer.start();
    // Advance simulation
    sim.advance(config.tmax);
    timer.stop();

    //Save expected values time series
    auto means   = sim.get_expected_values_time_series();
    auto moments = sim.get_moment_time_series(closure_order);

    int num_steps = static_cast<int>((config.tmax - init_time) / config.dt) + 1;
    std::vector<double> interpolation_tps(num_steps);
    for (int i = 0; i < num_steps; ++i) {
        interpolation_tps[i] = init_time + i * config.dt;
    }
    means     = mio::interpolate_simulation_result(means, interpolation_tps);
    auto done = means.export_csv(save_file + "means.csv");
    // Save moment time series
    auto moment_ts = mio::interpolate_simulation_result(moments.first, interpolation_tps);
    done           = moment_ts.export_csv(save_file + "moments.csv", moments.second);

    // Save total time
    Eigen::VectorXd time = Eigen::VectorXd::Constant(1, mio::timing::time_in_seconds(timer.get_elapsed_time()));
    mio::TimeSeries<double> total_time(1);
    total_time.add_time_point(0., time);
    auto finished_time = total_time.export_csv(save_file + "total_time.csv", {"Runtime"});

    return 0;
}
