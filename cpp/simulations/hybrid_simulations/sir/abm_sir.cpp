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

#include "abm/common_abm_loggers.h"
#include "abm/infection_state.h"
#include "abm/location_type.h"
#include "abm/mobility_rules.h"
#include "abm/model.h"
#include "abm/infectivity_functions.h"
#include "abm/parameters.h"
#include "abm/personal_rng.h"
#include "abm/time.h"
#include "abm/virus_variant.h"
#include "memilio/config.h"
#include "memilio/data/analyze_result.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/utils/parameter_distributions.h"
#include "memilio/utils/random_number_generator.h"
#include "setup.cpp"
#include <vector>

void run_abm_sim(int sim_num)
{
    mio::set_log_level(mio::LogLevel::warn);
    //Create model with one age group an constant transmission rate
    auto model = mio::abm::Model(1, &mio::abm::constant_infectivity,
                                 std::vector<mio::abm::Model::MobilityRuleType>{&mio::abm::get_buried});
    // Set transmission rate to lambda
    model.parameters.get<mio::abm::VirusShedFactor>() = mio::ParameterDistributionConstant(setup::lambda);
    // Set time from Exposed to Non-Symptomatic to an exponential distribution with rate gamma
    model.parameters.get<mio::abm::TimeExposedToNoSymptoms>() = mio::ParameterDistributionExponential(setup::gamma);
    // Set all time from NoSymptoms to Recovered to zero so an agent directly recovers after exposure (which we use here as Infected)
    model.parameters.get<mio::abm::TimeInfectedNoSymptomsToRecovered>() = mio::ParameterDistributionConstant(0.);
    // Set Symptoms per NoSymptoms to zero so all agents recover directly
    model.parameters.get<mio::abm::SymptomsPerInfectedNoSymptoms>() = 0.;
    // Set aerosol transmission rates to zero
    model.parameters.get<mio::abm::AerosolTransmissionRates>() = 0.;
    // Add only one location which is assigned to all persons
    auto home_id = model.add_location(mio::abm::LocationType::Home);
    // First the infected persons are added (here having infection state E which is used as I compartment)
    double num_initially_infected = 0;
    while (num_initially_infected < setup::I0) {
        auto pid     = model.add_person(home_id, mio::AgeGroup(0));
        auto& person = model.get_person(pid);
        person.set_assigned_location(mio::abm::LocationType::Home, home_id, model.get_id());
        auto p_rng = mio::abm::PersonalRandomNumberGenerator(person);
        person.add_new_infection(mio::abm::Infection(
            p_rng, mio::abm::VirusVariant::Wildtype, person.get_age(), model.parameters, mio::abm::TimePoint(setup::t0),
            &mio::abm::constant_infectivity, mio::abm::InfectionState::Exposed));
        num_initially_infected += 1;
    }
    // Then create susceptible persons
    int num_remaining_agents = int(setup::total_population - num_initially_infected);
    for (auto p = 0; p < num_remaining_agents; ++p) {
        auto pid     = model.add_person(home_id, mio::AgeGroup(0));
        auto& person = model.get_person(pid);
        person.set_assigned_location(mio::abm::LocationType::Home, home_id, model.get_id());
    }
    // Create simulation
    auto tp_max = mio::abm::TimePoint(setup::t0) + mio::abm::days(setup::tmax);
    auto sim    = mio::abm::Simulation(mio::abm::TimePoint(setup::t0), std::move(model));
    // History object aggregates the number of agents per infection state and time point
    mio::History<mio::abm::TimeSeriesWriter, mio::abm::LogInfectionState> history{
        Eigen::Index(mio::abm::InfectionState::Count)};

    // Run the simulation until tmax with the history object.
    sim.advance(tp_max, history);

    std::string output_file =
        "/home/bick_ju/Documents/TemporalHybridModel/sim_outputs/ABM/" + std::to_string(sim_num) + "_comps.csv";

    int num_steps = static_cast<int>(setup::tmax / setup::dt) + 1;

    std::vector<double> interpolation_tps(num_steps);

    for (int i = 0; i < num_steps; ++i) {
        interpolation_tps[i] = i * setup::dt;
    }
    auto result = mio::interpolate_simulation_result(std::get<0>(history.get_log()), interpolation_tps)
                      .export_csv(output_file, {"S", "E", "I_NS", "I_Sy", "I_Sev", "I_Crit", "R", "D"});
}

int main()
{
    int num_sims = 500;
    for (int sim = 0; sim < num_sims; ++sim) {
        run_abm_sim(sim);
    }
}
