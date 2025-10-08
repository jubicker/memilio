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
#include "memilio/data/analyze_result.h"
#include "ode_sir/infection_state.h"
#include "ode_sir/model.h"
#include "ode_sir/parameters.h"
#include "memilio/compartments/simulation.h"
#include "setup.cpp"

int main()
{

    mio::set_log_level(mio::LogLevel::debug);

    mio::log_info("Simulating SIR; t={} ... {} with dt = {}.", setup::t0, setup::tmax, setup::dt);

    mio::osir::Model<ScalarType> model(1);

    model.populations[{mio::AgeGroup(0), mio::osir::InfectionState::Infected}]  = setup::I0;
    model.populations[{mio::AgeGroup(0), mio::osir::InfectionState::Recovered}] = 0;
    model.populations[{mio::AgeGroup(0), mio::osir::InfectionState::Susceptible}] =
        setup::total_population - model.populations[{mio::AgeGroup(0), mio::osir::InfectionState::Infected}] -
        model.populations[{mio::AgeGroup(0), mio::osir::InfectionState::Recovered}];
    //Set contact frequency to total population
    mio::ContactMatrixGroup& contact_matrix =
        model.parameters.get<mio::osir::ContactPatterns<ScalarType>>().get_cont_freq_mat();
    contact_matrix[0].get_baseline().setConstant(setup::total_population);

    //rho = lambda if contact frequency is total population
    model.parameters.set<mio::osir::TransmissionProbabilityOnContact<ScalarType>>(setup::lambda);
    model.parameters.set<mio::osir::TimeInfected<ScalarType>>(1. / setup::gamma);

    model.check_constraints();

    auto sir = mio::simulate(setup::t0, setup::tmax, setup::dt, model);

    std::string output_file = "/home/bick_ju/Documents/TemporalHybridModel/sim_outputs/ODE/sir.csv";

    int num_steps = static_cast<int>(setup::tmax / setup::dt) + 1;

    std::vector<double> interpolation_tps(num_steps);

    for (int i = 0; i < num_steps; ++i) {
        interpolation_tps[i] = i * setup::dt;
    }
    auto result = mio::interpolate_simulation_result(sir, interpolation_tps).export_csv(output_file, {"S", "I", "R"});
}
