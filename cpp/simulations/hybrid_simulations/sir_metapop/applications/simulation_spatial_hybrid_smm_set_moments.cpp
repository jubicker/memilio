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
    // const size_t num_runs      = 10000;
    // double dt_switch           = 1.;
    const size_t closure_order = 3;
    const auto config          = Config::get_config(Config::ConfigType::Config2regionsk1);
    const size_t num_regions   = 2;
    // double min_step_size       = 0.0001;

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

    save_file += "/" + Config::switch_condition_string[1];
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

    // // Second temp hybrid model for second region and delete the rates per model in every region
    // mio::hybrid::SpatialHybridSimulation<num_regions, closure_order> spatial_hybrid_sim(config, num_runs, min_step_size,
    //                                                                                     dt_switch, config.t0);

    // mio::timing::BasicTimer timer;
    // timer.start();
    // spatial_hybrid_sim.advance(config.tmax);
    // timer.stop();

    // //Calculate outputs per region
    // for (size_t r = 0; r < num_regions; ++r) {
    //     auto region_sim = spatial_hybrid_sim.get_sim_by_region(r);
    //     // Calculate moments and expected values of sim set
    //     auto means_smm       = region_sim.get_model1().get_mean();
    //     auto moments_smm     = region_sim.get_model1().get_moments();
    //     auto done            = means_smm.export_csv(save_file + std::to_string(r) + "_means_smm.csv");
    //     done                 = moments_smm.export_csv(save_file + std::to_string(r) + "_moments_smm.csv",
    //                                                   region_sim.get_model1().get_moment_names());
    //     auto expected_values = region_sim.get_model2().get_expected_values_time_series();
    //     auto moments_moments = region_sim.get_model2().get_moment_time_series(closure_order);
    //     done = expected_values.export_csv(save_file + std::to_string(r) + "_expected_values_moments.csv");
    //     done = moments_moments.first.export_csv(save_file + std::to_string(r) + "_moments_moments.csv",
    //                                             moments_moments.second);

    //     // Get percentiles of smm simulation set
    //     std::vector<std::vector<mio::TimeSeries<double>>> sim_results;
    //     auto& all_results = region_sim.get_model1().get_result();
    //     //size_t run        = 0;
    //     for (auto& res : all_results) {
    //         //(void)res.export_csv(save_file + std::to_string(run) + "_comps.csv");
    //         sim_results.push_back({res});
    //         //run += 1;
    //     }
    //     // Save percentiles
    //     auto p05      = mio::ensemble_percentile(sim_results, 0.05);
    //     auto p25      = mio::ensemble_percentile(sim_results, 0.25);
    //     auto p50      = mio::ensemble_percentile(sim_results, 0.5);
    //     auto p75      = mio::ensemble_percentile(sim_results, 0.75);
    //     auto p95      = mio::ensemble_percentile(sim_results, 0.95);
    //     auto finished = p05[0].export_csv(save_file + std::to_string(r) + "_smm_p05.csv");
    //     finished      = p25[0].export_csv(save_file + std::to_string(r) + "_smm_p25.csv");
    //     finished      = p50[0].export_csv(save_file + std::to_string(r) + "_smm_p50.csv");
    //     finished      = p75[0].export_csv(save_file + std::to_string(r) + "_smm_p75.csv");

    //     // Save merged time series
    //     auto hybrid_result_means   = mio::merge_time_series(means_smm, expected_values).value();
    //     auto hybrid_result_moments = mio::merge_time_series(moments_smm, moments_moments.first).value();

    //     int num_steps = static_cast<int>(config.tmax / config.dt) + 1;
    //     std::vector<double> interpolation_tps(num_steps);

    //     for (int i = 0; i < num_steps; ++i) {
    //         interpolation_tps[i] = i * config.dt;
    //     }
    //     hybrid_result_means   = mio::interpolate_simulation_result(hybrid_result_means, interpolation_tps);
    //     hybrid_result_moments = mio::interpolate_simulation_result(hybrid_result_moments, interpolation_tps);
    //     done                  = hybrid_result_means.export_csv(save_file + std::to_string(r) + "_means.csv");
    //     done = hybrid_result_moments.export_csv(save_file + std::to_string(r) + "_moments.csv", moments_moments.second);
    // }

    // mio::TimeSeries<double> total_time(1);
    // Eigen::VectorXd time = Eigen::VectorXd::Constant(1, mio::timing::time_in_seconds(timer.get_elapsed_time()));
    // total_time.add_time_point(0., time);
    // auto finished_time = total_time.export_csv(save_file + "total_time.csv", {"Runtime"});

    //std::cout << "Spatial-hybrid Elapsed time: " << timer.get_elapsed_time() << std::endl << std::flush;

    return 0;
}
