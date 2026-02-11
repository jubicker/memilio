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
#include "hybrid/simulation_set.h"
#include "memilio/config.h"
#include "ode_sir/infection_state.h"
#include "simulations/hybrid_simulations/sir_metapop/config/config.cpp"
#include "simulations/hybrid_simulations/sir_metapop/library/moment_helper.h"
#include "simulations/hybrid_simulations/sir_metapop/library/smm_helper.h"
#include "memilio/data/analyze_result.h"
#include "memilio/timer/basic_timer.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/time_series.h"
#include "smm/simulation.h"
#include "smm_moments/simulation.h"
#include "smm_moments/closure_functions.h"
#include "models/hybrid/conversion_functions.cpp"
#include <cstddef>
#include <cstdint>
#include <omp.h>
#include <vector>

/**
 * @brief Run one simulation of the temporal-hybrid smm-ode-sir model. The comps.csv is saved as result.
 * @tparam NumRegions Number of regions
 * @param[in] sim_num Simulation number used as seed for the simulation.
 * @param[in] save_file File simulation output timeseries is saved to.
 * @param[in] config Simulation config i.e. parameters and populations used to initialize the model.
 * @param[in] rel_switch_value If this percentage of the total population is infected, it is switched from smm to ode.
 * @param[in, out] timer_init Timer used to time initialization.
 * @param[in, out] timer_sim Timer used to time simulation.
 * @return Simulation time series i.e. merged time series of smm and ode.
 */
template <size_t NumRegions>
mio::TimeSeries<double> run_hybrid_sim(size_t sim_num, std::string save_file, const Config::Config& config,
                                       double rel_switch_value, mio::timing::BasicTimer& timer_init,
                                       mio::timing::BasicTimer& timer_sim)
{
    timer_init.start();
    // Initialize smm
    auto smm_model = smm_helper::initialize_model<NumRegions>(config);
    // Set seed
    smm_model.get_rng().seed({static_cast<uint32_t>(sim_num)});

    // Initialize expected values
    Eigen::Array<double, Eigen::Dynamic, 1> expected_values_init(NumRegions *
                                                                 static_cast<size_t>(mio::osir::InfectionState::Count));
    expected_values_init.setZero();
    // Initialize moment expected values with smm populations
    for (size_t r = 0; r < NumRegions; ++r) {
        for (size_t s = 0; s < static_cast<size_t>(mio::osir::InfectionState::Count); ++s) {
            expected_values_init[r * static_cast<size_t>(mio::osir::InfectionState::Count) + s] =
                smm_model.populations[{mio::regions::Region(r), mio::osir::InfectionState(s)}];
        }
    }
    // Moment are all zero
    MomentArray<static_cast<size_t>(mio::osir::InfectionState::Count), NumRegions, 1> moments_array;
    moments_array.moments()[moments_array.flatten_index({0, 0, 0})] = 1.0;
    //Initialize moment model
    auto moment_model =
        moment_helper::initialize_model<NumRegions, 2>(expected_values_init, moments_array.moments(), config);

    // Create simulations
    auto sim_smm =
        mio::smm::Simulation<ScalarType, NumRegions, mio::osir::InfectionState>(smm_model, config.t0, config.dt);
    auto sim_moments = mio::smm_moments::Simulation<NumRegions, 2>(moment_model, config.t0, config.dt);

    // Define result functions
    const auto result_fct_smm = [](const mio::smm::Simulation<ScalarType, NumRegions, mio::osir::InfectionState>& sim,
                                   double /*t*/) {
        return sim.get_result();
    };
    const auto result_fct_moments = [](const mio::smm_moments::Simulation<NumRegions, 2>& sim, double /*t*/) {
        return sim.get_expected_values_time_series();
    };

    // Define switching conditions - Model switches when total number of infected is bigger that given percentage of the total population
    const auto condition = [rel_switch_value, &config](const mio::TimeSeries<double>& result_smm,
                                                       const mio::TimeSeries<double>& result_ode, bool smm_used) {
        double total_population = 0;
        for (size_t r = 0; r < NumRegions; ++r) {
            total_population += config.total_populations[r];
        }
        if (smm_used) {
            auto& last_value      = result_smm.get_last_value().eval();
            double total_infected = 0;
            for (size_t r = 0; r < NumRegions; ++r) {
                total_infected +=
                    last_value[r * (int)mio::osir::InfectionState::Count + (int)mio::osir::InfectionState::Infected];
            }
            if ((total_infected > rel_switch_value * total_population) || (total_infected < 1)) {
                return true;
            }
        }
        else {
            auto& last_value      = result_ode.get_last_value().eval();
            double total_infected = 0;
            for (size_t r = 0; r < NumRegions; ++r) {
                total_infected +=
                    last_value[r * (int)mio::osir::InfectionState::Count + (int)mio::osir::InfectionState::Infected];
            }
            if ((total_infected <= rel_switch_value * total_population) && (total_infected >= 1)) {
                return true;
            }
        }
        return false;
    };

    // Create hybrid simulation
    double dt_switch = 0.2;
    mio::hybrid::TemporalHybridSimulation<decltype(sim_smm), decltype(sim_moments), mio::TimeSeries<double>,
                                          mio::TimeSeries<double>>
        hybrid_sim(std::move(sim_smm), std::move(sim_moments), result_fct_smm, result_fct_moments, true, config.t0,
                   dt_switch);
    timer_init.stop();

    timer_sim.start();
    hybrid_sim.advance(config.tmax, condition);
    timer_sim.stop();

    // If both time series contain the same time point (which is the case at the switching time point and at simulation start), the values of the first model are taken
    auto hybrid_result = mio::merge_time_series(hybrid_sim.get_result_model1(), hybrid_sim.get_result_model2()).value();
    std::string output_file = save_file + std::to_string(sim_num) + "_comps.csv";

    // Interpolate to dt time steps
    int num_steps = static_cast<int>(config.tmax / config.dt) + 1;
    std::vector<double> interpolation_tps(num_steps);

    for (int i = 0; i < num_steps; ++i) {
        interpolation_tps[i] = i * config.dt;
    }
    auto result = mio::interpolate_simulation_result(hybrid_result, interpolation_tps);

    auto done = result.export_csv(output_file);

    return result;
}

int main()
{
    mio::set_log_level(mio::LogLevel::warn);
    const size_t num_runs         = 1;
    double dt_switch              = 1.;
    const size_t max_order        = 3;
    const auto config             = Config::get_config(Config::ConfigType::Config1);
    const size_t num_regions      = 1;
    const double rel_switch_value = 1.0;
    if (num_regions != config.num_regions) {
        mio::log_error("Number of regions doesn't match number of regions in config.");
    }

    std::string save_file = Config::SAVE_DIR + "Hybrid1/";
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

    auto smm_model = smm_helper::initialize_model<num_regions>(config);
    // Initialize expected values
    Eigen::Array<double, Eigen::Dynamic, 1> expected_values_init(num_regions *
                                                                 static_cast<size_t>(mio::osir::InfectionState::Count));
    expected_values_init.setZero();
    // Initialize moment expected values with smm populations
    for (size_t r = 0; r < num_regions; ++r) {
        for (size_t s = 0; s < static_cast<size_t>(mio::osir::InfectionState::Count); ++s) {
            expected_values_init[r * static_cast<size_t>(mio::osir::InfectionState::Count) + s] =
                smm_model.populations[{mio::regions::Region(r), mio::osir::InfectionState(s)}];
        }
    }
    // Moment are all zero
    MomentArray<static_cast<size_t>(mio::osir::InfectionState::Count), num_regions, 2> moments_array;
    moments_array.moments()[moments_array.flatten_index({0, 0, 0})] = 1.0;
    //Initialize moment model
    auto moment_model =
        moment_helper::initialize_model<num_regions, 2>(expected_values_init, moments_array.moments(), config);

    // Define result functions
    const auto result_fct_smm = [](const mio::smm::Simulation<ScalarType, num_regions, mio::osir::InfectionState>& sim,
                                   double /*t*/) {
        return sim.get_result();
    };
    const auto result_fct_moments = [](const mio::smm_moments::Simulation<num_regions, 2>& sim, double /*t*/) {
        return sim.get_expected_values_time_series();
    };

    // Define switching conditions - Model switches when total number of infected is bigger that given percentage of the total population
    const auto condition = [rel_switch_value, &config](const mio::TimeSeries<double>& result_smm,
                                                       const mio::TimeSeries<double>& /*result_ode*/, bool smm_used) {
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
        // else {
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

    auto sim_set = mio::hybrid::SimulationSet<num_regions, mio::osir::InfectionState, max_order>(
        num_runs, smm_model, moment_model, result_fct_smm, result_fct_moments, config.t0, config.dt, dt_switch);

    mio::timing::BasicTimer timer;
    timer.start();
    // Advance simulation set
    sim_set.advance(config.tmax, condition);
    timer.stop();

    // Convert result so they fit structure for ensemble_percentile fct
    std::vector<std::vector<mio::TimeSeries<double>>> sim_results;
    auto& all_results = sim_set.get_result();
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
    auto finished = p05[0].export_csv(save_file + "p05.csv");
    finished      = p25[0].export_csv(save_file + "p25.csv");
    finished      = p50[0].export_csv(save_file + "p50.csv");
    finished      = p75[0].export_csv(save_file + "p75.csv");
    finished      = p95[0].export_csv(save_file + "p95.csv");
    // Save means and moments
    finished = sim_set.get_mean().export_csv(save_file + "means.csv");
    finished = sim_set.get_moments().export_csv(save_file + "moments.csv", sim_set.get_moment_names());
    // Save times
    mio::TimeSeries<double> time_ts(1);
    auto& all_times = sim_set.get_sim_times();
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

    return 0;
}
