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
#include "memilio/timer/definitions.h"
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
    mio::set_log_level(mio::LogLevel::err);
    const size_t num_runs      = 1000;
    double dt_exchange         = 1.;
    const size_t closure_order = 3;
    auto config                = Config::get_config(Config::ConfigType::ConfigInfluenzaRegions_full);
    const size_t num_regions   = 4;
    double min_step_size       = 0.0001;

    if (num_regions != config.num_regions) {
        mio::log_error("Number of regions doesn't match number of regions in config.");
    }

    // INPUT - SIR/SIRS
    std::string save_file  = Config::SAVE_DIR + "Spatial-Hybrid2/";
    auto created_directory = mio::create_directory(save_file);
    save_file += config.name;

    created_directory = mio::create_directory(save_file);
    if (!created_directory) {
        printf("%s\n", created_directory.error().formatted_message().c_str());
        return -1;
    }

    save_file += "/combined_region";
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
    Condition.set_absolute_switch_threshold(100.);
    Condition.set_mean_gradient_threshold(0);
    Condition.set_R0_threshold(1.0);
    double first                     = 0.4;
    double first_bav                 = 0.3;
    double second                    = 0.95;
    double second_bav                = 0.95;
    double third                     = 0.8;
    double third_bav                 = 0.7;
    double fourth                    = 0.7;
    double fourth_bav                = 0.6;
    double fifth                     = 0.85;
    double fifth_bav                 = 0.85;
    double sixth                     = 0.85;
    double sixth_bav                 = 0.85;
    std::vector<double> damping_vec1 = {first, first, first_bav, first};
    std::vector<double> damping_vec2 = {second / first, second / first, second_bav / first_bav, second / first};
    std::vector<double> damping_vec3 = {third / second, third / second, third_bav / second_bav, third / second};
    std::vector<double> damping_vec4 = {fourth / third, fourth / third, fourth_bav / third_bav, fourth / third};
    std::vector<double> damping_vec5 = {fifth / fourth, fifth / fourth, fifth_bav / fourth_bav, fifth / fourth};
    std::vector<double> damping_vec6 = {sixth / fifth, sixth / fifth, sixth_bav / fifth_bav, sixth / fifth};

    mio::timing::BasicTimer timer;
    timer.start();
    spatial_hybrid_sim.advance(1323, Condition.combined_region, true);
    spatial_hybrid_sim.apply_dampings(damping_vec1);
    spatial_hybrid_sim.advance(1353, Condition.combined_region, true);
    spatial_hybrid_sim.apply_dampings(damping_vec2);
    spatial_hybrid_sim.advance(1554, Condition.combined_region, true);
    spatial_hybrid_sim.apply_dampings(damping_vec3);
    spatial_hybrid_sim.advance(1577, Condition.combined_region, true);
    spatial_hybrid_sim.apply_dampings(damping_vec4);
    spatial_hybrid_sim.advance(1680, Condition.combined_region, true);
    spatial_hybrid_sim.apply_dampings(damping_vec5);
    spatial_hybrid_sim.advance(1736, Condition.combined_region, true);
    spatial_hybrid_sim.apply_dampings(damping_vec6);
    spatial_hybrid_sim.advance(config.tmax, Condition.combined_region, true);
    timer.stop();

    // Save stochastic outputs
    auto done    = spatial_hybrid_sim.get_stochastic_mean(config.dt).export_csv(save_file + "stochastic_mean.csv");
    auto moments = spatial_hybrid_sim.get_stochastic_moments(config.dt);
    done         = moments.first.export_csv(save_file + "stochastic_moments.csv", moments.second);
    // Save deterministic outputs
    done    = spatial_hybrid_sim.get_deterministic_mean(config.dt).export_csv(save_file + "deterministic_mean.csv");
    moments = spatial_hybrid_sim.get_deterministic_moments(config.dt);
    done    = moments.first.export_csv(save_file + "deterministic_moments.csv", moments.second);
    // Save joint outputs
    done    = spatial_hybrid_sim.get_joint_mean(config.dt).export_csv(save_file + "means.csv");
    moments = spatial_hybrid_sim.get_joint_moments(config.dt);
    done    = moments.first.export_csv(save_file + "moments.csv", moments.second);
    // Save model used
    done = spatial_hybrid_sim.get_model_used_ts().export_csv(save_file + "model_used.csv");
    // Save runtime
    mio::TimeSeries<double> total_time(1);
    Eigen::VectorXd time = Eigen::VectorXd::Constant(1, mio::timing::time_in_seconds(timer.get_elapsed_time()));
    total_time.add_time_point(0., time);
    auto finished_time = total_time.export_csv(save_file + "total_time.csv", {"Runtime"});

    std::cout << "Spatial-hybrid Elapsed time: " << mio::timing::time_in_seconds(timer.get_elapsed_time()) << std::endl
              << std::flush;

    return 0;
}
