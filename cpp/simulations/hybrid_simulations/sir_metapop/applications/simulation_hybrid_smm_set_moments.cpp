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
#include "simulations/hybrid_simulations/sir_metapop/config/config.cpp"
#include "memilio/utils/logging.h"
#include "simulations/hybrid_simulations/sir_metapop/library/moment_helper.h"
#include "simulations/hybrid_simulations/sir_metapop/library/smm_helper.h"
#include "smm/simulation_set.h"
#include "smm_moments/simulation.h"
#include "models/hybrid/conversion_functions.cpp"

int main()
{
    mio::set_log_level(mio::LogLevel::warn);
    const size_t num_runs         = 10000;
    const size_t max_order        = 5;
    const auto config             = Config::get_config(Config::ConfigType::Config2);
    const size_t num_regions      = 1;
    const double rel_switch_value = 0.0;
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

    save_file += "/switch_value_" + std::to_string(rel_switch_value);
    created_directory = mio::create_directory(save_file);
    if (!created_directory) {
        printf("%s\n", created_directory.error().formatted_message().c_str());
        return -1;
    }
    save_file += "/";

    // Initialize smm
    auto smm_model = smm_helper::initialize_model<num_regions>(config);
    // Initialize smm simulation set
    auto sim_set = mio::smm::SimulationSet<num_regions, mio::osir::InfectionState, max_order>(num_runs, smm_model,
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
    MomentArray<static_cast<size_t>(mio::osir::InfectionState::Count), num_regions, max_order> moments_array;
    moments_array.moments()[moments_array.flatten_index({0, 0, 0})] = 1.0;
    auto moment_model =
        moment_helper::initialize_model<num_regions, max_order>(expected_values_init, moments_array.moments(), config);

    // Initialize moment simulation
    auto sim_moments = mio::smm_moments::Simulation<num_regions, max_order>(moment_model, config.t0, config.dt);
    // Set maximum dt of integrator to interpolation time points
    sim_moments.get_integrator_core().get_dt_max() = config.dt;

    // Define result functions
    const auto result_fct_smm = [](mio::smm::SimulationSet<num_regions, mio::osir::InfectionState, max_order>& sim,
                                   double /*t*/) {
        return sim.get_mean();
    };
    const auto result_fct_moments = [](mio::smm_moments::Simulation<num_regions, max_order>& sim, double /*t*/) {
        return sim.get_expected_values_time_series();
    };

    // Define switching conditions - Model switches when total number of infected is bigger that given percentage of the total population
    const auto condition = [rel_switch_value, &config](mio::TimeSeries<double>& result_smm,
                                                       mio::TimeSeries<double>& /*result_ode*/, bool smm_used) {
        double total_population = 0;
        for (size_t r = 0; r < num_regions; ++r) {
            total_population += config.total_populations[r];
        }
        if (smm_used) {
            auto& last_value      = result_smm.get_last_value().eval();
            double total_infected = 0;
            for (size_t r = 0; r < num_regions; ++r) {
                total_infected +=
                    last_value[r * (int)mio::osir::InfectionState::Count + (int)mio::osir::InfectionState::Infected];
            }
            if ((total_infected > rel_switch_value * total_population) || (total_infected < 1)) {
                return true;
            }
        }
        // else { //Switch back
        //     auto& last_value      = result_ode.get_last_value().eval();
        //     double total_infected = 0;
        //     for (size_t r = 0; r < num_regions; ++r) {
        //         total_infected +=
        //             last_value[r * (int)mio::osir::InfectionState::Count + (int)mio::osir::InfectionState::Infected];
        //     }
        //     if ((total_infected <= rel_switch_value * total_population) && (total_infected >= 1)) {
        //         return true;
        //     }
        // }
        return false;
    };

    double hybrid_check_dt = 0.5;
    mio::hybrid::TemporalHybridSimulation<decltype(sim_set), decltype(sim_moments), mio::TimeSeries<double>,
                                          mio::TimeSeries<double>>
        hybrid_sim(std::move(sim_set), std::move(sim_moments), result_fct_smm, result_fct_moments, true, config.t0,
                   hybrid_check_dt);

    mio::timing::BasicTimer timer;
    timer.start();
    hybrid_sim.advance(config.tmax, condition);
    timer.stop();
    std::cout << "Time: " << timer.get_elapsed_time() << std::endl;

    // Calculate moments and expected values of sim set
    auto means_smm   = hybrid_sim.get_model1().get_mean();
    auto moments_smm = hybrid_sim.get_model1().get_moments();
    auto done        = means_smm.export_csv(save_file + "means_smm.csv");
    done = moments_smm.export_csv(save_file + "moments_smm.csv", hybrid_sim.get_model1().get_moment_names());
    auto expected_values = hybrid_sim.get_model2().get_expected_values_time_series();
    auto moments_moments = hybrid_sim.get_model2().get_moment_time_series(max_order);
    done                 = expected_values.export_csv(save_file + "expected_values_moments.csv");
    done                 = moments_moments.first.export_csv(save_file + "moments_moments.csv", moments_moments.second);

    auto hybrid_result_means   = mio::merge_time_series(means_smm, expected_values).value();
    auto hybrid_result_moments = mio::merge_time_series(moments_smm, moments_moments.first).value();

    int num_steps = static_cast<int>(config.tmax / config.dt) + 1;
    std::vector<double> interpolation_tps(num_steps);

    for (int i = 0; i < num_steps; ++i) {
        interpolation_tps[i] = i * config.dt;
    }
    hybrid_result_means   = mio::interpolate_simulation_result(hybrid_result_means, interpolation_tps);
    hybrid_result_moments = mio::interpolate_simulation_result(hybrid_result_moments, interpolation_tps);
    done = hybrid_result_means.export_csv(save_file + "merged_expected_values.csv", {"muS_r0", "muI_r0", "muR_r0"});
    done = hybrid_result_moments.export_csv(save_file + "merged_moments.csv", moments_moments.second);

    return 0;
}
