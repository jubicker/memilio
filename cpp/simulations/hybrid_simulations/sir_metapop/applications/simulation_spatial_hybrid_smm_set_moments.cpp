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

#include "hybrid/spatial_hybrid_model.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/uncertain_value.h"
#include "ode_sir/infection_state.h"
#include "simulations/hybrid_simulations/sir_metapop/config/config.cpp"
#include "memilio/utils/logging.h"
#include "simulations/hybrid_simulations/sir_metapop/library/condition_functions.h"
#include "simulations/hybrid_simulations/sir_metapop/library/moment_helper.h"
#include "simulations/hybrid_simulations/sir_metapop/library/smm_helper.h"
#include "smm/simulation_set.h"
#include "smm_moments/simulation.h"
#include "smm_moments/closure_functions.h"
#include "smm/model.h"
#include <cstddef>
#include <limits>
#include <utility>

int main()
{
    mio::set_log_level(mio::LogLevel::warn);
    const size_t num_runs      = 10000;
    double dt_exchange         = 1.;
    const size_t closure_order = 3;
    const auto config          = Config::get_config(Config::ConfigType::Config3);
    const size_t num_regions   = 1;
    double min_step_size       = 0.0001;

    if (num_regions != config.num_regions) {
        mio::log_error("Number of regions doesn't match number of regions in config.");
    }

    std::string save_file = Config::SAVE_DIR + "Spatial-Hybrid2/";
    save_file += config.name;

    auto created_directory = mio::create_directory(save_file);
    if (!created_directory) {
        printf("%s\n", created_directory.error().formatted_message().c_str());
        return -1;
    }

    save_file += "/pure_stochastic";
    created_directory = mio::create_directory(save_file);
    if (!created_directory) {
        printf("%s\n", created_directory.error().formatted_message().c_str());
        return -1;
    }

    save_file += "/" + Config::closure_string[0];
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

    save_file += "/";

    // Create spatial-hybrid simulation - INPUT: Closure function
    mio::hybrid::SpatialHybridSimulation<num_regions, closure_order> spatial_hybrid_sim(
        config, num_runs, min_step_size, dt_exchange,
        &mio::smm_moments::truncation_closure<num_regions, closure_order>);

    // Create switching condition - INPUT: threshold value etc. for condition
    mio::hybrid::SwitchingCondition<num_regions, closure_order> Condition;
    Condition.set_config(config);
    Condition.set_mean_stddev_relation(0.3);
    Condition.set_var_gradient_threshold(-1.);

    mio::timing::BasicTimer timer;
    timer.start();
    spatial_hybrid_sim.advance(config.tmax, Condition.pure_stochastic, true);
    timer.stop();

    // Save stochastic outputs
    double result_dt = 0.1;
    auto done        = spatial_hybrid_sim.get_stochastic_mean(result_dt).export_csv(save_file + "stochastic_mean.csv");
    auto moments     = spatial_hybrid_sim.get_stochastic_moments(result_dt);
    done             = moments.first.export_csv(save_file + "stochastic_moments.csv", moments.second);
    // Save deterministic outputs
    done    = spatial_hybrid_sim.get_deterministic_mean(result_dt).export_csv(save_file + "deterministic_mean.csv");
    moments = spatial_hybrid_sim.get_deterministic_moments(result_dt);
    done    = moments.first.export_csv(save_file + "deterministic_moments.csv", moments.second);
    // Save joint outputs
    done    = spatial_hybrid_sim.get_joint_mean(result_dt).export_csv(save_file + "joint_mean.csv");
    moments = spatial_hybrid_sim.get_joint_moments(result_dt);
    done    = moments.first.export_csv(save_file + "joint_moments.csv", moments.second);
    // Save model used
    done = spatial_hybrid_sim.get_model_used_ts().export_csv(save_file + "model_used.csv");
    // Save runtime
    mio::TimeSeries<double> total_time(1);
    Eigen::VectorXd time = Eigen::VectorXd::Constant(1, mio::timing::time_in_seconds(timer.get_elapsed_time()));
    total_time.add_time_point(0., time);
    auto finished_time = total_time.export_csv(save_file + "total_time.csv", {"Runtime"});

    std::cout << "Spatial-hybrid Elapsed time: " << timer.get_elapsed_time() << std::endl << std::flush;

    return 0;
}
