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

#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include "ode_sir/infection_state.h"
#include "smm/simulation.h"
#include "smm/model.h"
#include "smm/parameters.h"
#include "memilio/data/analyze_result.h"
#include "memilio/epidemiology/adoption_rate.h"
#include "hybrid_simulations/sir/moment_equation_array.h"
#include "setup.cpp"
#include <boost/math/special_functions/math_fwd.hpp>
#include <string>
#include <vector>

mio::TimeSeries<double> run_smm_sim(int sim_num, std::string save_file)
{
    using Model = mio::smm::Model<1, mio::osir::InfectionState>;

    Model model;
    model.populations[{mio::regions::Region(0), mio::osir::InfectionState::Infected}]  = setup::I0;
    model.populations[{mio::regions::Region(0), mio::osir::InfectionState::Recovered}] = 0;
    model.populations[{mio::regions::Region(0), mio::osir::InfectionState::Susceptible}] =
        setup::total_population - model.populations[{mio::regions::Region(0), mio::osir::InfectionState::Infected}] -
        model.populations[{mio::regions::Region(0), mio::osir::InfectionState::Recovered}];

    // Second-order adoption rate is lambda * N
    std::vector<mio::AdoptionRate<mio::osir::InfectionState>> adoption_rates;
    // lambda has to be multiplied with the contact frequency which is the total population
    adoption_rates.push_back({mio::osir::InfectionState::Susceptible,
                              mio::osir::InfectionState::Infected,
                              mio::regions::Region(0),
                              setup::lambda * setup::total_population,
                              {{mio::osir::InfectionState::Infected, 1.}}});
    adoption_rates.push_back({mio::osir::InfectionState::Infected,
                              mio::osir::InfectionState::Recovered,
                              mio::regions::Region(0),
                              setup::gamma,
                              {}});

    model.parameters.get<mio::smm::AdoptionRates<mio::osir::InfectionState>>() = adoption_rates;

    auto sim = mio::smm::Simulation(model, setup::t0, setup::dt);
    sim.advance(setup::tmax);

    std::string output_file = save_file + std::to_string(sim_num) + "_comps.csv";

    int num_steps = static_cast<int>(setup::tmax / setup::dt) + 1;

    std::vector<double> interpolation_tps(num_steps);

    for (int i = 0; i < num_steps; ++i) {
        interpolation_tps[i] = i * setup::dt;
    }
    auto result = mio::interpolate_simulation_result(sim.get_result(), interpolation_tps);
    auto done   = result.export_csv(output_file, {"S", "I", "R"});
    return result;
}

int main()
{
    int num_sims               = 1000;
    const size_t closing_order = 3;
    SIR::Moments<closing_order> moments({setup::total_population - setup::I0, setup::I0, 0});
    std::vector<mio::TimeSeries<double>> sim_results;
    std::string save_file = "/home/bick_ju/Documents/TemporalHybridModel/sim_outputs/SMM/I0=" +
                            std::to_string(static_cast<int>(setup::I0)) + "/";
    setup::save_setup(save_file + "setup.csv");

    // Run multiple simulations
    for (int sim = 0; sim < num_sims; ++sim) {
        sim_results.push_back(run_smm_sim(sim, save_file));
    }

    // Calculate moments from the simulation results
    mio::TimeSeries<ScalarType> moments_ts = moments.calculate_moments_from_sim(sim_results);
    auto result_moments                    = moments.get_moment_ts_and_names(moments_ts, closing_order);
    auto finished = result_moments.first.export_csv(save_file + "moments.csv", result_moments.second);
}
