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

#include "hybrid/temporal_hybrid_model.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/uncertain_value.h"
#include "simulations/hybrid_simulations/sir_metapop/config/config.cpp"
#include "memilio/utils/logging.h"
#include "simulations/hybrid_simulations/sir_metapop/library/moment_helper.h"
#include "simulations/hybrid_simulations/sir_metapop/library/smm_helper.h"
#include "simulations/hybrid_simulations/sir_metapop/library/condition_functions.h"
#include "smm/simulation_set.h"
#include "smm_moments/simulation.h"
#include "models/hybrid/conversion_functions.cpp"
#include "smm_moments/closure_functions.h"
#include <cstddef>
#include <limits>

int main()
{
    mio::set_log_level(mio::LogLevel::warn);
    const size_t num_runs         = 10000;
    double dt_switch              = 1.;
    const size_t closure_order    = 3;
    const auto config             = Config::sir::get_config(Config::sir::ConfigType::Config5);
    const size_t num_regions      = 1;
    const double rel_switch_value = 0.3;
    double min_step_size          = 0.0001;
    size_t closure                = 0;
    size_t condition              = 3;
    auto closure_func             = &mio::smm_moments::truncation_closure<num_regions, closure_order>;
    if (closure == 0) {
        closure_func = &mio::smm_moments::truncation_closure<num_regions, closure_order>;
    }
    else if (closure == 1) {
        closure_func = &mio::smm_moments::pairapprox_closure<num_regions, closure_order>;
    }
    else if (closure == 2) {
        closure_func = &mio::smm_moments::lognormal_closure<num_regions, closure_order>;
    }
    else {
        mio::log_error("Unkown closure type: ", closure);
    }

    if (num_regions != config.num_regions) {
        mio::log_error("Number of regions doesn't match number of regions in config.");
    }

    std::string save_file = Config::SAVE_DIR + "Hybrid2/";
    save_file += config.name;

    auto created_directory = mio::create_directory(save_file);
    if (!created_directory) {
        printf("%s\n", created_directory.error().formatted_message().c_str());
        return -1;
    }

    save_file += "/" + Config::switch_condition_string[condition];
    created_directory = mio::create_directory(save_file);
    if (!created_directory) {
        printf("%s\n", created_directory.error().formatted_message().c_str());
        return -1;
    }

    if (condition == 0) {
        save_file += "/" + std::to_string(rel_switch_value);
        created_directory = mio::create_directory(save_file);
        if (!created_directory) {
            printf("%s\n", created_directory.error().formatted_message().c_str());
            return -1;
        }
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

    save_file += "/";

    // Initialize smm
    auto smm_model = smm_helper::initialize_model<num_regions>(config);
    // Initialize smm simulation set
    auto sim_set = mio::smm::SimulationSet<num_regions, mio::osir::InfectionState, closure_order>(num_runs, smm_model,
                                                                                                  config.t0, config.dt);
    // Initialize moment model
    // Initial expected values
    Eigen::Array<double, Eigen::Dynamic, 1> expected_values_init(num_regions *
                                                                 static_cast<size_t>(mio::osir::InfectionState::Count));
    expected_values_init.setZero();
    for (size_t r = 0; r < num_regions; ++r) {
        for (size_t s = 0; s < static_cast<size_t>(mio::osir::InfectionState::Count); ++s) {
            expected_values_init[r * static_cast<size_t>(mio::osir::InfectionState::Count) + s] =
                smm_model.populations[{mio::regions::Region(r), mio::osir::InfectionState(s)}];
        }
    }
    // Initial moments are all zero
    MomentArray<static_cast<size_t>(mio::osir::InfectionState::Count), num_regions, closure_order> moments_array;
    std::array<int, static_cast<size_t>(mio::osir::InfectionState::Count) * num_regions> zero_index;
    zero_index.fill(0);
    moments_array.moments()[moments_array.flatten_index(zero_index)] = 1.0;
    auto moment_model = moment_helper::initialize_model<num_regions, closure_order>(
        expected_values_init, moments_array.moments(), config, closure_func);

    // Initialize moment simulation
    auto sim_moments = mio::smm_moments::Simulation<num_regions, closure_order>(moment_model, config.t0, config.dt);
    // Set maximum dt of integrator to interpolation time points
    sim_moments.get_integrator_core().get_dt_max() = config.dt;
    if (min_step_size > 0) {
        sim_moments.get_integrator_core().get_dt_min() = min_step_size;
    }

    // Set switching condition - INPUT
    mio::hybrid::SwitchingCondition<num_regions, closure_order, decltype(config)> Condition;
    Condition.set_config(config);
    //Condition.set_rel_switch_threshold(rel_switch_value);
    Condition.set_mean_stddev_relation(0.3);

    // Initialize hybrid simulation with correct result functions for switching condition - INPUT
    mio::hybrid::TemporalHybridSimulation<decltype(sim_set), decltype(sim_moments), std::vector<double>,
                                          std::vector<double>>
        hybrid_sim(std::move(sim_set), std::move(sim_moments), Condition.current_smm_relations,
                   Condition.current_moment_relations, true, config.t0, dt_switch);

    mio::timing::BasicTimer timer;
    timer.start();
    // Run hybrid simulation with switching condition - INPUT
    hybrid_sim.advance(config.tmax, Condition.mean_stddev_relation_condition);
    timer.stop();

    // Calculate moments and expected values of sim set
    auto means_smm   = hybrid_sim.get_model1().get_mean();
    auto moments_smm = hybrid_sim.get_model1().get_moments();
    auto done        = means_smm.export_csv(save_file + "means_smm.csv");
    done = moments_smm.export_csv(save_file + "moments_smm.csv", hybrid_sim.get_model1().get_moment_names());
    auto expected_values = hybrid_sim.get_model2().get_expected_values_time_series();
    auto moments_moments = hybrid_sim.get_model2().get_moment_time_series(closure_order);
    done                 = expected_values.export_csv(save_file + "expected_values_moments.csv");
    done                 = moments_moments.first.export_csv(save_file + "moments_moments.csv", moments_moments.second);

    // Get percentiles of smm simulation set
    std::vector<std::vector<mio::TimeSeries<double>>> sim_results;
    auto& all_results = hybrid_sim.get_model1().get_result();
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
    auto finished = p05[0].export_csv(save_file + "smm_p05.csv");
    finished      = p25[0].export_csv(save_file + "smm_p25.csv");
    finished      = p50[0].export_csv(save_file + "smm_p50.csv");
    finished      = p75[0].export_csv(save_file + "smm_p75.csv");
    finished      = p95[0].export_csv(save_file + "smm_p95.csv");

    // Save merged time series
    auto hybrid_result_means   = mio::merge_time_series(means_smm, expected_values, true).value();
    auto hybrid_result_moments = mio::merge_time_series(moments_smm, moments_moments.first, true).value();

    int num_steps = static_cast<int>(config.tmax / config.dt) + 1;
    std::vector<double> interpolation_tps(num_steps);

    for (int i = 0; i < num_steps; ++i) {
        interpolation_tps[i] = i * config.dt;
    }
    hybrid_result_means   = mio::interpolate_simulation_result(hybrid_result_means, interpolation_tps);
    hybrid_result_moments = mio::interpolate_simulation_result(hybrid_result_moments, interpolation_tps);
    done                  = hybrid_result_means.export_csv(save_file + "means.csv");
    done                  = hybrid_result_moments.export_csv(save_file + "moments.csv", moments_moments.second);

    // Save times
    mio::TimeSeries<double> time_ts(1);
    auto& all_times = hybrid_sim.get_model1().get_sim_times();
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

    std::cout << "Tempral-hybrid Elapsed time: " << timer.get_elapsed_time() << std::endl << std::flush;

    return 0;
}
