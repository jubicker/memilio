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

#include "hybrid/conversion_functions.h"
#include "d_abm/single_well.h"
#include "memilio/geography/regions.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/random_number_generator.h"
#include "ode_secir/infection_state.h"
#include "memilio/epidemiology/age_group.h"
#include "simulations/hybrid_simulations/sir_metapop/library/moment_helper.h"
#include <algorithm>
#include <cstddef>
#include <numeric>

namespace mio
{
namespace hybrid
{
template <>
void convert_model(dabm::Simulation<SingleWell<mio::osecir::InfectionState>>& current_model,
                   smm::Simulation<ScalarType, 1, mio::osecir::InfectionState>& target_model)
{
    auto& current_result = current_model.get_result();
    auto& target_result  = target_model.get_result();
    if (current_result.get_last_time() < target_result.get_last_time()) {
        mio::log_error(
            "Conversion from dabm to smm not possible because last smm time point is bigger than last dabm time point");
    }
    if (target_result.get_last_time() < current_result.get_last_time()) {
        target_result.add_time_point(current_result.get_last_time());
    }
    // Update result timeseries
    target_result.get_last_value() = current_result.get_last_value();
    // Update model populations
    for (int i = 0; i < (int)mio::osecir::InfectionState::Count; ++i) {
        target_model.get_model().populations[{regions::Region(0), mio::osecir::InfectionState(i)}] =
            current_result.get_last_value()[target_model.get_model().populations.get_flat_index(
                {regions::Region(0), mio::osecir::InfectionState(i)})];
    }
}

template <>
void convert_model(smm::Simulation<ScalarType, 1, mio::osecir::InfectionState>& current_model,
                   dabm::Simulation<SingleWell<mio::osecir::InfectionState>>& target_model)
{
    auto& current_result = current_model.get_result();
    auto& target_result  = target_model.get_result();
    if (current_result.get_last_time() < target_result.get_last_time()) {
        mio::log_error("Conversion from smm to dabm not possible because last dabm time point is bigger than last smm "
                       "time point.");
    }
    if (target_result.get_last_time() < current_result.get_last_time()) {
        target_result.add_time_point(current_result.get_last_time());
    }
    // Update result timeseries
    target_result.get_last_value() = current_result.get_last_value();

    // Update agents' infection state and sample agents position
    auto current_pop = current_result.get_last_value().eval();
    double total_pop = std::accumulate(current_pop.begin(), current_pop.end(), 0.0);
    SWPositionSampler pos_rng{{-1, -1}, {1, 1}, 0.1};
    auto& state_rng = DiscreteDistribution<size_t>::get_instance();
    auto& abm_pop   = target_model.get_model().populations;
    if (abm_pop.size() == 0) {
        log_info("Diffusive ABM does not contain any agents. Population is initialized from SMM population.");
        int num_agents = 0;
        while (num_agents < total_pop) {
            auto position        = pos_rng();
            auto infection_state = state_rng(thread_local_rng(), current_pop);
            abm_pop.push_back(
                SingleWell<mio::osecir::InfectionState>::Agent{position, mio::osecir::InfectionState(infection_state)});
            current_pop[infection_state] = std::max(0., current_pop[infection_state] - 1);
            num_agents++;
        }
    }
    else {
        assert(int(total_pop) == int(abm_pop.size()) && "Population sizes of dabm and smm do not match.");
        for (auto& a : abm_pop) {
            auto infection_state         = state_rng(thread_local_rng(), current_pop);
            a.position                   = pos_rng();
            a.status                     = mio::osecir::InfectionState(infection_state);
            current_pop[infection_state] = std::max(0., current_pop[infection_state] - 1);
        }
    }
}

template <>
void convert_model(dabm::Simulation<SingleWell<mio::osecir::InfectionState>>& current_model,
                   mio::Simulation<ScalarType, mio::osecir::Model<ScalarType>>& target_model)
{
    auto& current_result = current_model.get_result();
    auto& target_result  = target_model.get_result();
    if (current_result.get_last_time() < target_result.get_last_time()) {
        mio::log_error("Conversion from dabm to ODE-SECIR not possible because last ODE-SECIR time point is bigger "
                       "than last dabm time point.");
    }
    if (target_result.get_last_time() < current_result.get_last_time()) {
        target_result.add_time_point(current_result.get_last_time());
    }
    // If the secir model has more than one age group, the compartments from the dabm are equally distributed to all age groups
    size_t num_age_groups = target_result.get_last_value().size() / (int)mio::osecir::InfectionState::Count;
    for (int i = 0; i < (int)mio::osecir::InfectionState::Count; ++i) {
        for (size_t age_group = 0; age_group < num_age_groups; ++age_group) {
            double pop_value = current_result.get_last_value()[(int)mio::osecir::InfectionState(i)] / num_age_groups;
            // Update result timeseries
            target_result.get_last_value()[target_model.get_model().populations.get_flat_index(
                {mio::AgeGroup(age_group), mio::osecir::InfectionState(i)})] = pop_value;
            // Update model populations
            target_model.get_model().populations[{mio::AgeGroup(age_group), mio::osecir::InfectionState(i)}] =
                pop_value;
        }
    }
}

template <>
void convert_model(mio::Simulation<ScalarType, mio::osecir::Model<ScalarType>>& current_model,
                   dabm::Simulation<SingleWell<mio::osecir::InfectionState>>& target_model)
{
    auto& current_result = current_model.get_result();
    auto& target_result  = target_model.get_result();
    if (current_result.get_last_time() < target_result.get_last_time()) {
        mio::log_error("Conversion from ODE-SECIR to dabm not possible because last dabm time point is bigger than "
                       "last ODE-SECIR time point.");
    }
    if (target_result.get_last_time() < current_result.get_last_time()) {
        target_result.add_time_point(current_result.get_last_time());
    }

    size_t num_age_groups = current_result.get_last_value().size() / (int)mio::osecir::InfectionState::Count;
    target_result.get_last_value().setZero();
    // Update dabm time series
    // ODE-SECIR model's age groups are aggregated as dabm does not have age groups
    for (size_t age_group = 0; age_group < num_age_groups; ++age_group) {
        for (int i = 0; i < (int)mio::osecir::InfectionState::Count; ++i) {
            target_result.get_last_value()[i] +=
                current_result.get_last_value()[current_model.get_model().populations.get_flat_index(
                    {mio::AgeGroup(age_group), mio::osecir::InfectionState(i)})];
        }
    }

    // Update agents' infection state and sample agents position
    auto current_pop = target_result.get_last_value().eval();
    double total_pop = std::accumulate(current_pop.begin(), current_pop.end(), 0.0);
    SWPositionSampler pos_rng{{-1, -1}, {1, 1}, 0.1};
    auto& state_rng = DiscreteDistribution<size_t>::get_instance();
    auto& abm_pop   = target_model.get_model().populations;
    if (abm_pop.size() == 0) {
        log_info("Diffusive ABM does not contain any agents. Population is initialized from ODE-SECIR compartments.");
        int num_agents = 0;
        while (num_agents < total_pop) {
            auto position        = pos_rng();
            auto infection_state = state_rng(thread_local_rng(), current_pop);
            abm_pop.push_back(
                SingleWell<mio::osecir::InfectionState>::Agent{position, mio::osecir::InfectionState(infection_state)});
            current_pop[infection_state] = std::max(0., current_pop[infection_state] - 1);
            num_agents++;
        }
    }
    else {
        assert(int(total_pop) == int(abm_pop.size()) && "Population sizes of dabm and ODE-SECIR do not match.");
        for (auto& a : abm_pop) {
            auto infection_state         = state_rng(thread_local_rng(), current_pop);
            a.position                   = pos_rng();
            a.status                     = mio::osecir::InfectionState(infection_state);
            current_pop[infection_state] = std::max(0., current_pop[infection_state] - 1);
        }
    }
}

template <>
void convert_model(smm::Simulation<double, 1, mio::osir::InfectionState>& current_model,
                   smm_moments::Simulation<1, 2>& target_model)
{
    auto& current_result = current_model.get_result();
    auto& target_result  = target_model.get_result();
    if (current_result.get_last_time() < target_result.get_last_time()) {
        mio::log_error("Conversion from smm to dabm not possible because last dabm time point is bigger than last smm "
                       "time point.");
    }
    if (target_result.get_last_time() < current_result.get_last_time()) {
        target_result.add_time_point(current_result.get_last_time());
    }

    // Update result timeseries
    auto smm_values     = current_result.get_last_value();
    auto moments_values = target_result.get_last_value();
    for (auto i = 0; i < smm_values.size(); ++i) {
        // Set expected values
        moments_values[i] = smm_values[i];
    }

    // Update model populations
    for (int i = 0; i < (int)mio::osir::InfectionState::Count; ++i) {
        target_model.get_model().populations[{regions::Region(0), mio::osir::InfectionState(i)}] =
            current_result.get_last_value()[target_model.get_model().populations.get_flat_index(
                {regions::Region(0), mio::osir::InfectionState(i)})];
    }
}

template <>
void convert_model(smm_moments::Simulation<1, 2>& current_model,
                   smm::Simulation<double, 1, mio::osir::InfectionState>& target_model)
{
    auto& current_result = current_model.get_result();
    auto& target_result  = target_model.get_result();
    if (current_result.get_last_time() < target_result.get_last_time()) {
        mio::log_error("Conversion from smm to dabm not possible because last dabm time point is bigger than last smm "
                       "time point.");
    }
    if (target_result.get_last_time() < current_result.get_last_time()) {
        target_result.add_time_point(current_result.get_last_time());
    }

    // Update result timeseries
    auto moments_values = current_result.get_last_value();
    auto smm_values     = target_result.get_last_value();
    for (auto i = 0; i < smm_values.size(); ++i) {
        // Set expected values
        smm_values[i] = moments_values[i];
    }

    // Update model populations
    for (int i = 0; i < (int)mio::osir::InfectionState::Count; ++i) {
        target_model.get_model().populations[{regions::Region(0), mio::osir::InfectionState(i)}] =
            current_result.get_last_value()[target_model.get_model().populations.get_flat_index(
                {regions::Region(0), mio::osir::InfectionState(i)})];
    }
}

template <>
void convert_model(smm::Simulation<double, 2, mio::osir::InfectionState>& current_model,
                   smm_moments::Simulation<2, 2>& target_model)
{
    auto& current_result = current_model.get_result();
    auto& target_result  = target_model.get_result();
    if (current_result.get_last_time() < target_result.get_last_time()) {
        mio::log_error(
            "Conversion from smm to moments not possible because last moment time point is bigger than last smm "
            "time point.");
    }
    if (target_result.get_last_time() < current_result.get_last_time()) {
        target_result.add_time_point(current_result.get_last_time());
    }

    // Update result timeseries
    auto smm_values     = current_result.get_last_value();
    auto moments_values = target_result.get_last_value();
    for (auto i = 0; i < smm_values.size(); ++i) {
        // Set expected values
        moments_values[i] = smm_values[i];
    }

    // Update model populations
    for (int i = 0; i < (int)mio::osir::InfectionState::Count; ++i) {
        target_model.get_model().populations[{regions::Region(0), mio::osir::InfectionState(i)}] =
            current_result.get_last_value()[target_model.get_model().populations.get_flat_index(
                {regions::Region(0), mio::osir::InfectionState(i)})];
        target_model.get_model().populations[{regions::Region(1), mio::osir::InfectionState(i)}] =
            current_result.get_last_value()[target_model.get_model().populations.get_flat_index(
                {regions::Region(0), mio::osir::InfectionState(i)})];
    }
}

template <>
void convert_model(smm_moments::Simulation<2, 2>& current_model,
                   smm::Simulation<double, 2, mio::osir::InfectionState>& target_model)
{
    auto& current_result = current_model.get_result();
    auto& target_result  = target_model.get_result();
    if (current_result.get_last_time() < target_result.get_last_time()) {
        mio::log_error(
            "Conversion from smm to moment not possible because last moment time point is bigger than last smm "
            "time point.");
    }
    if (target_result.get_last_time() < current_result.get_last_time()) {
        target_result.add_time_point(current_result.get_last_time());
    }

    // Update result timeseries
    auto moments_values = current_result.get_last_value();
    auto smm_values     = target_result.get_last_value();
    for (auto i = 0; i < smm_values.size(); ++i) {
        // Set expected values
        smm_values[i] = moments_values[i];
    }

    // Update model populations
    for (int i = 0; i < (int)mio::osir::InfectionState::Count; ++i) {
        target_model.get_model().populations[{regions::Region(0), mio::osir::InfectionState(i)}] =
            current_result.get_last_value()[target_model.get_model().populations.get_flat_index(
                {regions::Region(0), mio::osir::InfectionState(i)})];
        target_model.get_model().populations[{regions::Region(0), mio::osir::InfectionState(i)}] =
            current_result.get_last_value()[target_model.get_model().populations.get_flat_index(
                {regions::Region(1), mio::osir::InfectionState(i)})];
    }
}

/**
 * @brief Specialization for SMM Set -> Moment model.
 * Conversion does the following steps:
 * - Copy last mean and moment values from smm set to moment result.
 * - Set last value in smm simulation results (attribute in smm::SimulationSet and the individual ones in smm::Simulation ) to zero.
 * - Set all smm populations to zero.
 * - Set last values of smm set means and moments to zero.
 */
template <>
void convert_model(smm::SimulationSet<1, mio::osir::InfectionState, 3>& current_model,
                   smm_moments::Simulation<1, 3>& target_model)
{
    // Calculate smm statistics
    auto& smm_means     = current_model.get_mean();
    auto& smm_moments   = current_model.get_moments();
    auto moment_names   = current_model.get_moment_names();
    auto& target_result = target_model.get_result();

    if (smm_means.get_last_time() < target_result.get_last_time()) {
        mio::log_error("Conversion from smm simulation set to moments not possible because last moment time point is "
                       "bigger than last smm sim set time point.");
    }

    if (target_result.get_last_time() < smm_means.get_last_time()) {
        target_result.add_time_point(smm_means.get_last_time());
    }

    auto last_moment_values = target_result.get_last_value();

    // Copy means to moment model
    auto last_smm_mean_values = smm_means.get_last_value();
    for (auto i = 0; i < last_smm_mean_values.size(); ++i) {
        last_moment_values[i] = last_smm_mean_values[i];
    }
    // Copy moments to moments model
    auto last_smm_moment_values = smm_moments.get_last_value();
    for (size_t i = 0; i < moment_names.size(); ++i) {
        const auto moment_index =
            moment_helper::moment_to_indices<static_cast<size_t>(mio::osir::InfectionState::Count), 1>(moment_names[i]);
        size_t flat_index = target_model.get_model().moments.flatten_index(moment_index) +
                            target_model.get_model().populations.get_num_compartments();
        last_moment_values[flat_index] = last_smm_moment_values[i];
    }

    // Set last value in smm simulation results to zero
    auto& smm_set_result = current_model.get_result();
    auto& smm_simulation = current_model.get_simulations();
    for (size_t sim = 0; sim < smm_simulation.size(); ++sim) {
        smm_simulation[sim].get_result().get_last_value().setZero();
        smm_set_result[sim].get_last_value().setZero();
        smm_simulation[sim].get_model().populations.array().setZero();
    }

    // Set last means and moments in smm simulation set to zero
    smm_means.get_last_value().setZero();
    smm_moments.get_last_value().setZero();
}

/**
 * @brief Specialization for Moment model -> SMM Set.
 * Conversion does the following steps:
 * - Sample number of agents for every SMM simulation (normally distributed with mean and variance from moment model)
 * - Set populations in SMM model to sampled values.
 * - Set last value in smm simulation results (attribute in smm::SimulationSet and the individual ones in smm::Simulation ) to sampled values.
 * - Recalculate last means and moments in SMM set.
 * - Set last means and moments in moment sim to 0.
 */
template <>
void convert_model(smm_moments::Simulation<1, 3>& current_model,
                   smm::SimulationSet<1, mio::osir::InfectionState, 3>& target_model)
{
    auto mean_ts        = current_model.get_expected_values_time_series();
    auto means          = mean_ts.get_last_value();
    auto moment_results = current_model.get_result().get_last_value();

    // Set mean and moments
    if (mean_ts.get_last_time() < target_model.get_mean().get_last_time()) {
        mio::log_error("Conversion from moments to smm simulation set not possible because last smm simulation set "
                       "time point is bigger than last moment time point.");
    }

    auto& smm_set_sims = target_model.get_simulations();
    auto& smm_set_res  = target_model.get_result();
    for (size_t sim = 0; sim < smm_set_sims.size(); ++sim) {
        auto& sim_result = smm_set_sims[sim].get_result();
        auto& sim_pop    = smm_set_sims[sim].get_model().populations;
        for (size_t comp = 0; comp < static_cast<size_t>(mio::osir::InfectionState::Count); ++comp) {
            double mean = means[comp];
            std::array<int, static_cast<size_t>(mio::osir::InfectionState::Count)> indices_var;
            indices_var.fill(0);
            indices_var[comp] = 2;
            double var        = moment_results[current_model.get_model().populations.get_num_compartments() +
                                        current_model.get_model().moments.flatten_index(indices_var)];
            // Sample number of agents for simulation
            double sim_value = std::max(0., std::round(mio::NormalDistribution<double>::get_instance()(
                                                smm_set_sims[sim].get_model().get_rng(), mean, var)));
            // Set population in Simulation object
            sim_pop[{mio::regions::Region(0), mio::osir::InfectionState(comp)}] = sim_value;

            if (sim_result.get_last_time() < mean_ts.get_last_time()) {
                sim_result.add_time_point(mean_ts.get_last_time());
            }
            // Set last result in Simulation object
            sim_result.get_last_value()[comp] = sim_value;

            if (smm_set_res[sim].get_last_time() < mean_ts.get_last_time()) {
                smm_set_res[sim].add_time_point(mean_ts.get_last_time());
            }
            // Set last result in SimulationSet object
            smm_set_res[sim].get_last_value()[comp] = sim_value;
        }
    }

    // Recalculate last means and moments in smm simulation set
    target_model.recalculate_last_mean();
    target_model.recalculate_last_moments();

    // Set last means and moments in moment sim to 0
    moment_results.setZero();
}

// Two regions, closure order 3
template <>
void convert_model(smm::SimulationSet<2, mio::osir::InfectionState, 3>& current_model,
                   smm_moments::Simulation<2, 3>& target_model)
{
    // Calculate smm statistics
    auto& smm_means     = current_model.get_mean();
    auto& smm_moments   = current_model.get_moments();
    auto moment_names   = current_model.get_moment_names();
    auto& target_result = target_model.get_result();

    if (smm_means.get_last_time() < target_result.get_last_time()) {
        mio::log_error("Conversion from smm simulation set to moments not possible because last moment time point is "
                       "bigger than last smm sim set time point.");
    }

    if (target_result.get_last_time() < smm_means.get_last_time()) {
        target_result.add_time_point(smm_means.get_last_time());
    }

    auto last_moment_values = target_result.get_last_value();

    // Copy means to moment model
    auto last_smm_mean_values = smm_means.get_last_value();
    for (auto i = 0; i < last_smm_mean_values.size(); ++i) {
        last_moment_values[i] = last_smm_mean_values[i];
    }
    // Copy moments to moments model
    auto last_smm_moment_values = smm_moments.get_last_value();
    for (size_t i = 0; i < moment_names.size(); ++i) {
        const auto moment_index =
            moment_helper::moment_to_indices<static_cast<size_t>(mio::osir::InfectionState::Count), 2>(moment_names[i]);
        size_t flat_index = target_model.get_model().moments.flatten_index(moment_index) +
                            target_model.get_model().populations.get_num_compartments();
        last_moment_values[flat_index] = last_smm_moment_values[i];
    }
}

//TODO
template <>
void convert_model(smm_moments::Simulation<2, 3>& current_model,
                   smm::SimulationSet<2, mio::osir::InfectionState, 3>& target_model)
{
    const size_t num_regions = 2;
    auto mean_ts             = current_model.get_expected_values_time_series();
    auto means               = mean_ts.get_last_value();
    auto moment_ts           = current_model.get_moment_time_series(target_model.get_order()).first;
    // Set mean and moments
    if (mean_ts.get_last_time() < target_model.get_mean().get_last_time()) {
        mio::log_error("Conversion from moments to smm simulation set not possible because last smm simulation set "
                       "time point is bigger than last moment time point.");
    }

    if (target_model.get_mean().get_last_time() < mean_ts.get_last_time()) {
        target_model.get_mean().add_time_point(mean_ts.get_last_time());
        target_model.get_moments().add_time_point(moment_ts.get_last_time());
    }

    target_model.get_mean().get_last_value()    = means;
    target_model.get_moments().get_last_value() = moment_ts.get_last_value();

    // Set current state of all simulations of target model to mean
    auto& target_simulations = target_model.get_simulations();
    auto& target_results     = target_model.get_result();
    for (size_t sim = 0; sim < target_simulations.size(); ++sim) {
        if (mean_ts.get_last_time() < target_simulations[sim].get_result().get_last_time()) {
            mio::log_error("Conversion from moments to smm simulation set not possible because last smm simulation set "
                           "time point is bigger than last moment time point.");
        }

        if (target_simulations[sim].get_result().get_last_time() < mean_ts.get_last_time()) {
            target_simulations[sim].get_result().add_time_point(mean_ts.get_last_time());
            target_results[sim].add_time_point(mean_ts.get_last_time());
        }
        auto smm_values = target_simulations[sim].get_result().get_last_value();
        // Set expected values
        smm_values                           = means;
        target_results[sim].get_last_value() = means;

        // Update smm populations
        for (size_t region = 0; region < num_regions; ++region) {
            for (int i = 0; i < (int)mio::osir::InfectionState::Count; ++i) {
                target_simulations[sim]
                    .get_model()
                    .populations[{regions::Region(region), mio::osir::InfectionState(i)}] =
                    means[target_simulations[sim].get_model().populations.get_flat_index(
                        {regions::Region(region), mio::osir::InfectionState(i)})];
            }
        }
    }
    target_model.advance_t(mean_ts.get_last_time());
    target_model.advance_t_index();
}

// One region, closure order 5
template <>
void convert_model(smm::SimulationSet<1, mio::osir::InfectionState, 5>& current_model,
                   smm_moments::Simulation<1, 5>& target_model)
{
    // Calculate smm statistics
    auto& smm_means     = current_model.get_mean();
    auto& smm_moments   = current_model.get_moments();
    auto moment_names   = current_model.get_moment_names();
    auto& target_result = target_model.get_result();

    if (smm_means.get_last_time() < target_result.get_last_time()) {
        mio::log_error("Conversion from smm simulation set to moments not possible because last moment time point is "
                       "bigger than last smm sim set time point.");
    }

    if (target_result.get_last_time() < smm_means.get_last_time()) {
        target_result.add_time_point(smm_means.get_last_time());
    }

    auto last_moment_values = target_result.get_last_value();

    // Copy means to moment model
    auto last_smm_mean_values = smm_means.get_last_value();
    for (auto i = 0; i < last_smm_mean_values.size(); ++i) {
        last_moment_values[i] = last_smm_mean_values[i];
    }
    // Copy moments to moments model
    auto last_smm_moment_values = smm_moments.get_last_value();
    for (size_t i = 0; i < moment_names.size(); ++i) {
        const auto moment_index =
            moment_helper::moment_to_indices<static_cast<size_t>(mio::osir::InfectionState::Count), 1>(moment_names[i]);
        size_t flat_index = target_model.get_model().moments.flatten_index(moment_index) +
                            target_model.get_model().populations.get_num_compartments();
        last_moment_values[flat_index] = last_smm_moment_values[i];
    }
}

template <>
void convert_model(smm_moments::Simulation<1, 5>& current_model,
                   smm::SimulationSet<1, mio::osir::InfectionState, 5>& target_model)
{
    auto mean_ts   = current_model.get_expected_values_time_series();
    auto means     = mean_ts.get_last_value();
    auto moment_ts = current_model.get_moment_time_series(target_model.get_order()).first;
    // Set mean and moments
    if (mean_ts.get_last_time() < target_model.get_mean().get_last_time()) {
        mio::log_error("Conversion from moments to smm simulation set not possible because last smm simulation set "
                       "time point is bigger than last moment time point.");
    }

    if (target_model.get_mean().get_last_time() < mean_ts.get_last_time()) {
        target_model.get_mean().add_time_point(mean_ts.get_last_time());
        target_model.get_moments().add_time_point(moment_ts.get_last_time());
    }

    target_model.get_mean().get_last_value()    = means;
    target_model.get_moments().get_last_value() = moment_ts.get_last_value();

    // Set current state of all simulations of target model to mean
    auto& target_simulations = target_model.get_simulations();
    auto& target_results     = target_model.get_result();
    for (size_t sim = 0; sim < target_simulations.size(); ++sim) {
        if (mean_ts.get_last_time() < target_simulations[sim].get_result().get_last_time()) {
            mio::log_error("Conversion from moments to smm simulation set not possible because last smm simulation set "
                           "time point is bigger than last moment time point.");
        }

        if (target_simulations[sim].get_result().get_last_time() < mean_ts.get_last_time()) {
            target_simulations[sim].get_result().add_time_point(mean_ts.get_last_time());
            target_results[sim].add_time_point(mean_ts.get_last_time());
        }
        auto smm_values = target_simulations[sim].get_result().get_last_value();
        // Set expected values
        smm_values                           = means;
        target_results[sim].get_last_value() = means;

        // Update smm populations
        for (int i = 0; i < (int)mio::osir::InfectionState::Count; ++i) {
            target_simulations[sim].get_model().populations[{regions::Region(0), mio::osir::InfectionState(i)}] =
                means[target_simulations[sim].get_model().populations.get_flat_index(
                    {regions::Region(0), mio::osir::InfectionState(i)})];
        }
    }
    target_model.advance_t(mean_ts.get_last_time());
    target_model.advance_t_index();
}

} //namespace hybrid

} //namespace mio
