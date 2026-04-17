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

#ifndef MIO_SPATIAL_HYBRID_MODEL_H
#define MIO_SPATIAL_HYBRID_MODEL_H

#include "hybrid/temporal_hybrid_model.h"
#include "memilio/utils/random_number_generator.h"
#include "ode_sir/infection_state.h"
#include "simulations/hybrid_simulations/sir_metapop/library/moment_helper.h"
#include "simulations/hybrid_simulations/sir_metapop/library/moment_helper.h"
#include "smm/simulation_set.h"
#include "smm_moments/parameters.h"
#include "smm_moments/simulation.h"
#include "smm_moments/closure_functions.h"
#include "smm/model.h"
#include "simulations/hybrid_simulations/sir_metapop/config/config.h"
#include "simulations/hybrid_simulations/sir_metapop/library/smm_helper.h"
#include "hybrid/conversion_functions.h"
#include "hybrid/exchange_functions.h"
#include <cassert>
#include <cmath>
#include <cstddef>
#include <memory>
#include <numeric>
#include <vector>
#include <unordered_set>

namespace mio
{
namespace hybrid
{

template <size_t num_regions, size_t closure_order>
class SpatialHybridSimulation
{

public:
    using SMMSetSim           = mio::smm::SimulationSet<num_regions, mio::osir::InfectionState, closure_order>;
    using MomentSim           = mio::smm_moments::Simulation<num_regions, closure_order>;
    using ClosureFunctionType = typename mio::smm_moments::Model<num_regions, closure_order>::ClosureFunctionType;
    using switching_condition =
        std::function<bool(SMMSetSim& current_stochastic_model, MomentSim& current_deterministic_model,
                           bool stochastic_used, size_t region)>;

    SpatialHybridSimulation(
        const Config::Config& config, size_t num_runs, double min_step_size, double dt,
        ClosureFunctionType closure_func = &mio::smm_moments::truncation_closure<num_regions, closure_order>)
        : m_stochastic_simulation(initialize_stochastic_model(num_runs))
        , m_moment_simulation(initialize_deterministic_model(closure_func, min_step_size))
        , m_stochastic_regions(num_regions)
        , m_deterministic_regions(0)
        , m_t(config.t0)
        , m_dt(dt)
        , m_config(std::make_shared<Config::Config>(config))
    {
        // Initially all regions are modeled stochastically
        std::iota(m_stochastic_regions.begin(), m_stochastic_regions.end(), 0);
    }

    /**
     * @brief Advance spatial-hybrid simulation.
     * @param[in] tmax Simulation end time point.
     * @param[in] condition Switching condition for determining when to switch a region between stochastic and deterministic modeling.
     */
    void advance(double tmax, const switching_condition& condition)
    {
        while (m_t < tmax) {
            auto& stochastic_model_means   = m_stochastic_simulation.get_mean();
            auto& stochastic_model_moments = m_stochastic_simulation.get_moments();
            auto& moment_model_result      = m_moment_simulation.get_result();
            assert(m_stochastic_simulation.get_mean().get_last_time() == moment_model_result.get_last_time());
            // Evaluate switching conditions for all stochastic regions
            std::unordered_set<size_t> regions_to_switch_to_deterministic;
            for (size_t region : m_stochastic_regions) {
                if (condition(m_stochastic_simulation, m_moment_simulation, true, region)) {
                    regions_to_switch_to_deterministic.insert(region);
                }
            }
            if (regions_to_switch_to_deterministic.size() >
                0) { // Process regions that switch from stochastic to deterministic
                // Add regions that switched to deterministic to moment simulation
                add_regions_to_moment_model(regions_to_switch_to_deterministic);
                // Set transition and adoption rates for switched regions in moment model
                set_rates_in_moment_model(regions_to_switch_to_deterministic);
                // Remove switched regions from stochastic model
                remove_regions_from_stochastic_model(regions_to_switch_to_deterministic);
                // Delete transition and adoption rates for switched regions in stochastic model
                remove_rates_from_stochastic_model(regions_to_switch_to_deterministic);
            }
            // Evaluate switching conditions for all deterministic regions
            std::unordered_set<size_t> regions_to_switch_to_stochastic;
            for (size_t region : m_deterministic_regions) {
                if (condition(m_stochastic_simulation, m_moment_simulation, false, region)) {
                    regions_to_switch_to_stochastic.insert(region);
                }
            }
            if (regions_to_switch_to_stochastic.size() >
                0) { // Process regions that switch from deterministic to stochastic
                // Add regions that switched to stochastic to stochastic simulation
                add_regions_to_stochastic_model(regions_to_switch_to_stochastic);
                // Set transition and adoption rates for switched regions in stochastic model
                set_rates_in_stochastic_model(regions_to_switch_to_stochastic);
                // Remove switched regions from moment model
                remove_regions_from_moment_model(regions_to_switch_to_stochastic);
                // Delete transition and adoption rates for switched regions in moment model
                remove_rates_from_moment_model(regions_to_switch_to_stochastic);
            }
            { // Adapt vectors for stochastic and deterministic regions
                // Add new deterministic regions to moment regions vector and remove from stochastic regions vector
                m_deterministic_regions.insert(m_deterministic_regions.end(),
                                               regions_to_switch_to_deterministic.begin(),
                                               regions_to_switch_to_deterministic.end());
                m_stochastic_regions.erase(std::remove_if(m_stochastic_regions.begin(), m_stochastic_regions.end(),
                                                          [&](size_t r) {
                                                              return regions_to_switch_to_deterministic.count(r) > 0;
                                                          }),
                                           m_stochastic_regions.end());
                // Add new stochastic regions to stochastic regions vector and remove from deterministic regions vector
                m_stochastic_regions.insert(m_stochastic_regions.end(), regions_to_switch_to_stochastic.begin(),
                                            regions_to_switch_to_stochastic.end());
                m_deterministic_regions.erase(std::remove_if(m_deterministic_regions.begin(),
                                                             m_deterministic_regions.end(),
                                                             [&](size_t r) {
                                                                 return regions_to_switch_to_stochastic.count(r) > 0;
                                                             }),
                                              m_deterministic_regions.end());
            }
            { // Advance simulations for all regions
                double next_step = std::min(tmax, m_t + m_dt);
                // Advance stochastic simulation
                m_stochastic_simulation.advance(next_step);
                // Advance moment simulation
                m_moment_simulation.advance(next_step);
                m_t = next_step;
            }
            { // Exchange agents
                if (m_deterministic_regions.size() > 0 && m_stochastic_regions.size() > 0) {
                    exchange_moment_to_stochastic();
                }
                if (m_deterministic_regions.size() > 0 && m_stochastic_regions.size() > 0) {
                    exchange_stochastic_to_moment();
                }
            }
        }
    }

private:
    /**
     * @brief Check whether region is active in the given moment, meaning that it has entries > 0 in the multi-index.
     */
    bool is_region_active(
        const std::array<int, static_cast<size_t>(mio::osir::InfectionState::Count) * num_regions>& multiidx,
        size_t region)
    {
        const size_t base = region * static_cast<size_t>(mio::osir::InfectionState::Count);

        for (size_t s = 0; s < static_cast<size_t>(mio::osir::InfectionState::Count); ++s) {
            if (multiidx[base + s] > 0) {
                return true;
            }
        }
        return false;
    }

    /**
     * @brief Check whether moment has entries > 0 only in the allowed regions.
     */
    bool moment_only_in_regions(
        const std::array<int, static_cast<size_t>(mio::osir::InfectionState::Count) * num_regions>& multiidx,
        const std::unordered_set<size_t>& allowed_regions)
    {
        for (size_t r = 0; r < num_regions; ++r) {
            // Iterate over all regions and check whether the considered moment has entries > 0 for the current region
            if (is_region_active(multiidx, r)) {
                // Check if entries > 0 allowed in the current region
                if (allowed_regions.find(r) == allowed_regions.end()) {
                    return false;
                }
            }
        }
        return true;
    }

    void add_regions_to_moment_model(const std::unordered_set<size_t>& regions_to_switch)
    {
        if (regions_to_switch.empty()) {
            return;
        }
        auto& moment_model_result      = m_moment_simulation.get_result();
        auto& stochastic_model_means   = m_stochastic_simulation.get_mean();
        auto& stochastic_model_moments = m_stochastic_simulation.get_moments();
        for (size_t region : regions_to_switch) {
            // Copy means of switching regions to moment model
            for (size_t state = 0; state < static_cast<size_t>(mio::osir::InfectionState::Count); ++state) {
                assert(moment_model_result.get_last_value()[m_moment_simulation.get_model().populations.get_flat_index(
                           {mio::regions::Region(region), mio::osir::InfectionState(state)})] ==
                       0.); // Check that region is not already modeled in moment model
                moment_model_result.get_last_value()[m_moment_simulation.get_model().populations.get_flat_index(
                    {mio::regions::Region(region), mio::osir::InfectionState(state)})] =
                    stochastic_model_means
                        .get_last_value()[region * static_cast<size_t>(mio::osir::InfectionState::Count) + state];
            }
        }

        // Copy moments of switching regions to moment model
        for (size_t i = 0; i < m_stochastic_simulation.get_moment_indices().size(); ++i) {
            const auto moment_index = m_stochastic_simulation.get_moment_indices()[i];

            // Copy only moments that are partially or fully in switching regions
            bool switch_moment = moment_only_in_regions(moment_index, regions_to_switch);

            if (switch_moment) { // If all regions that have entries > 0 are switching, so we can switch the moment
                size_t flat_index = m_moment_simulation.get_model().moments.flatten_index(moment_index) +
                                    m_moment_simulation.get_model().populations.get_num_compartments();
                assert(moment_model_result.get_last_value()[flat_index] ==
                       0.); // Check that region is not already modeled in moment model
                moment_model_result.get_last_value()[flat_index] = stochastic_model_moments.get_last_value()[i];
            }
        }
    }

    void remove_regions_from_stochastic_model(const std::unordered_set<size_t>& regions_to_switch)
    {
        if (regions_to_switch.empty()) {
            return;
        }
        auto& stochastic_sims    = m_stochastic_simulation.get_simulations();
        auto& stochastic_results = m_stochastic_simulation.get_result();
        // Reset results and simulations of simulation set
        for (size_t sim = 0; sim < stochastic_sims.size(); ++sim) {
            auto& sim_result = stochastic_sims[sim].get_result();
            auto& sim_pop    = stochastic_sims[sim].get_model().populations;
            for (size_t region_to : regions_to_switch) {
                for (size_t comp = 0; comp < static_cast<size_t>(mio::osir::InfectionState::Count); ++comp) {
                    // Reset last result in all Simulation objects
                    sim_result
                        .get_last_value()[static_cast<size_t>(mio::osir::InfectionState::Count) * region_to + comp] = 0;
                    // Reset population in all Simulation objects
                    sim_pop[{mio::regions::Region(region_to), mio::osir::InfectionState(comp)}] = 0;
                    // Reset last result in SimulationSet object
                    stochastic_results[sim]
                        .get_last_value()[static_cast<size_t>(mio::osir::InfectionState::Count) * region_to + comp] = 0;
                }
            }
        }

        //Reset and update mean and moment ts of simulation set
        auto& stochastic_model_means   = m_stochastic_simulation.get_mean();
        auto& stochastic_model_moments = m_stochastic_simulation.get_moments();
        stochastic_model_means.remove_last_time_point();
        m_stochastic_simulation.recalculate_last_mean();
        stochastic_model_moments.remove_last_time_point();
        m_stochastic_simulation.recalculate_last_moments();
    }

    void set_rates_in_moment_model(const std::unordered_set<size_t>& regions_to_switch)
    {
        if (regions_to_switch.empty()) {
            return;
        }
        for (size_t region : regions_to_switch) {
            m_moment_simulation.get_model()
                .parameters.template get<mio::smm_moments::TransmissionRate>()[mio::regions::Region(region)] =
                m_config->lambdas[region];
            m_moment_simulation.get_model()
                .parameters.template get<mio::smm_moments::RecoveryRate>()[mio::regions::Region(region)] =
                m_config->gamma;

            for (auto& rate : m_config->transition_rates) {
                if (rate.from == mio::regions::Region(region)) {
                    m_moment_simulation.get_model().parameters.template get<mio::smm_moments::TransitionRate>()[{
                        rate.status, rate.from, rate.to}] = rate.factor;
                }
            }
        }
    }

    void remove_rates_from_stochastic_model(const std::unordered_set<size_t>& regions_to_switch)
    {
        if (regions_to_switch.empty()) {
            return;
        }
        for (size_t region : regions_to_switch) {
            m_stochastic_simulation.get_model()
                .parameters.template get<mio::smm_moments::TransmissionRate>()[mio::regions::Region(region)] = 0.0;
            m_stochastic_simulation.get_model()
                .parameters.template get<mio::smm_moments::RecoveryRate>()[mio::regions::Region(region)] = 0.0;

            for (auto& rate : m_config->transition_rates) {
                if (rate.from == mio::regions::Region(region)) {
                    m_stochastic_simulation.get_model().parameters.template get<mio::smm_moments::TransitionRate>()[{
                        rate.status, rate.from, rate.to}] = 0.0;
                }
            }
        }
    }

    void add_regions_to_stochastic_model(const std::unordered_set<size_t>& regions_to_switch)
    {
        if (regions_to_switch.empty()) {
            return;
        }
        auto moment_results      = m_moment_simulation.get_result().get_last_value();
        auto& stochastic_sims    = m_stochastic_simulation.get_simulations();
        auto& stochastic_results = m_stochastic_simulation.get_result();
        for (size_t region : regions_to_switch) {
            for (size_t state = 0; state < static_cast<size_t>(mio::osir::InfectionState::Count); ++state) {
                assert(m_stochastic_simulation.get_mean()
                           .get_last_value()[region * static_cast<size_t>(mio::osir::InfectionState::Count) + state] ==
                       0.); // Check that region is not already modeled in stochastic model
                for (size_t sim = 0; sim < stochastic_sims.size(); ++sim) {
                    auto& sim_result = stochastic_sims[sim].get_result();
                    auto& sim_pop    = stochastic_sims[sim].get_model().populations;

                    // Get current mean and variance from moment model for region and state
                    double mean =
                        moment_results[region * static_cast<size_t>(mio::osir::InfectionState::Count) + state];
                    std::array<int, static_cast<size_t>(mio::osir::InfectionState::Count) * num_regions> indices_var;
                    indices_var.fill(0);
                    indices_var[region * static_cast<size_t>(mio::osir::InfectionState::Count) + state] = 2;
                    double var = moment_results[m_moment_simulation.get_model().populations.get_num_compartments() +
                                                m_moment_simulation.get_model().moments.flatten_index(indices_var)];
                    // Sample number of agents for stochastic trajectory normally distributed
                    double sim_value =
                        std::max(0., std::round(mio::NormalDistribution<double>::get_instance()(
                                         stochastic_sims[sim].get_model().get_rng(), mean, std::sqrt(var))));
                    // Set population in Simulation object
                    sim_pop[{mio::regions::Region(region), mio::osir::InfectionState(state)}] = sim_value;

                    // Set sampled value in SMMSet result object
                    assert(mean_ts.get_last_time() == stochastic_results[sim].get_last_time());
                    stochastic_results
                        .get_last_value()[region * static_cast<size_t>(mio::osir::InfectionState::Count) + state] =
                        sim_value;
                    // Set sampled value in Simulation object
                    assert(mean_ts.get_last_time() == sim_result.get_last_time());
                    sim_result
                        .get_last_value()[region * static_cast<size_t>(mio::osir::InfectionState::Count) + state] =
                        sim_value;
                }
            }
        }

        //Reset and update mean and moment ts of simulation set
        auto& stochastic_model_means   = m_stochastic_simulation.get_mean();
        auto& stochastic_model_moments = m_stochastic_simulation.get_moments();
        stochastic_model_means.remove_last_time_point();
        m_stochastic_simulation.recalculate_last_mean();
        stochastic_model_moments.remove_last_time_point();
        m_stochastic_simulation.recalculate_last_moments();
    }

    void set_rates_in_stochastic_model(const std::unordered_set<size_t>& regions_to_switch)
    {
        if (regions_to_switch.empty()) {
            return;
        }
        for (size_t region : regions_to_switch) {
            m_stochastic_simulation.get_model()
                .parameters.template get<mio::smm_moments::TransmissionRate>()[mio::regions::Region(region)] =
                m_config->lambdas[region];
            m_stochastic_simulation.get_model()
                .parameters.template get<mio::smm_moments::RecoveryRate>()[mio::regions::Region(region)] =
                m_config->gamma;

            for (auto& rate : m_config->transition_rates) {
                if (rate.from == mio::regions::Region(region)) {
                    m_stochastic_simulation.get_model().parameters.template get<mio::smm_moments::TransitionRate>()[{
                        rate.status, rate.from, rate.to}] = rate.factor;
                }
            }
        }
    }

    void remove_regions_from_moment_model(const std::unordered_set<size_t>& regions_to_switch)
    {
        if (regions_to_switch.empty()) {
            return;
        }
        auto moment_result = m_moment_simulation.get_result();
        for (size_t region : regions_to_switch) {
            for (size_t index = 0; index < static_cast<size_t>(moment_result.get_num_elements()); ++index) {
                // Mean entries for corresponding regions are set to zero
                if (index >= region * static_cast<size_t>(mio::osir::InfectionState::Count) &&
                    index < region * static_cast<size_t>(mio::osir::InfectionState::Count) +
                                static_cast<size_t>(mio::osir::InfectionState::Count)) {
                    moment_result.get_last_value()[index] = 0;
                }
                if (index >= m_moment_simulation.get_model().populations.get_num_compartments()) {
                    auto multiindex = m_moment_simulation.get_model().moments.unflatten_index(
                        index - m_moment_simulation.get_model().populations.get_num_compartments());
                    int order_region = std::accumulate(
                        multiindex.begin() + static_cast<size_t>(mio::osir::InfectionState::Count) * region,
                        multiindex.begin() + static_cast<size_t>(mio::osir::InfectionState::Count) * region +
                            static_cast<size_t>(mio::osir::InfectionState::Count),
                        0.0);
                    if (order_region > 0) { // Moments with entries > 0 in switched region are set to zero
                        moment_result.get_last_value()[index] = 0;
                    }
                }
            }
        }
    }

    void remove_rates_from_moment_model(const std::unordered_set<size_t>& regions_to_switch)
    {
        if (regions_to_switch.empty()) {
            return;
        }
        for (size_t region : regions_to_switch) {
            m_moment_simulation.get_model()
                .parameters.template get<mio::smm_moments::TransmissionRate>()[mio::regions::Region(region)] = 0.0;
            m_moment_simulation.get_model()
                .parameters.template get<mio::smm_moments::RecoveryRate>()[mio::regions::Region(region)] = 0.0;

            for (auto& rate : m_config->transition_rates) {
                if (rate.from == mio::regions::Region(region)) {
                    m_moment_simulation.get_model().parameters.template get<mio::smm_moments::TransitionRate>()[{
                        rate.status, rate.from, rate.to}] = 0.0;
                }
            }
        }
    }

    void exchange_moment_to_stochastic()
    {
        auto moment_results      = m_moment_simulation.get_result();
        auto& stochastic_sims    = m_stochastic_simulation.get_simulations();
        auto& stochastic_results = m_stochastic_simulation.get_result();
        for (size_t region : m_stochastic_regions) {
            bool region_exchanged = false;
            for (size_t state = 0; state < static_cast<size_t>(mio::osir::InfectionState::Count); ++state) {
                double mean =
                    moment_results
                        .get_last_value()[region * static_cast<size_t>(mio::osir::InfectionState::Count) + state];
                if (mean < 1.) { // No (or not enough) agents to exchange for this region and state
                    continue;
                }
                else {
                    region_exchanged = true;
                    std::array<int, static_cast<size_t>(mio::osir::InfectionState::Count) * num_regions> indices_var;
                    indices_var.fill(0);
                    indices_var[region * static_cast<size_t>(mio::osir::InfectionState::Count) + state] = 2;
                    double var =
                        moment_results
                            .get_last_value()[m_moment_simulation.get_model().populations.get_num_compartments() +
                                              m_moment_simulation.get_model().moments.flatten_index(indices_var)];
                    for (size_t sim = 0; sim < stochastic_sims.size(); ++sim) {
                        auto& sim_result = stochastic_sims[sim].get_result();
                        auto& sim_pop    = stochastic_sims[sim].get_model().populations;
                        // Sample number of incoming agents for simulation
                        double sim_value =
                            std::max(0., std::round(mio::NormalDistribution<double>::get_instance()(
                                             stochastic_sims[sim].get_model().get_rng(), mean, std::sqrt(var))));

                        // Assert time last time point of moment and stochastic simulation is the same
                        assert(m_moment_simulation.get_result().get_last_time() == sim_result.get_last_time());
                        assert(m_moment_simulation.get_result().get_last_time() ==
                               stochastic_results[sim].get_last_time());
                        // Add agents to last result in all Simulation objects
                        sim_result
                            .get_last_value()[region * static_cast<size_t>(mio::osir::InfectionState::Count) + state] +=
                            sim_value;
                        // Add agents to population in all Simulation objects
                        sim_pop[{mio::regions::Region(region), mio::osir::InfectionState(state)}] += sim_value;
                        // Add agents to last result in SimulationSet object
                        stochastic_results[sim]
                            .get_last_value()[static_cast<size_t>(mio::osir::InfectionState::Count) * region + state] +=
                            sim_value;
                    }
                }
            }

            if (region_exchanged) {
                // If region was exchanged, set means and moments in stochastic regions in moment model to zero
                for (size_t index = 0; index < static_cast<size_t>(moment_results.get_num_elements()); ++index) {
                    // Means
                    if (index >= region * static_cast<size_t>(mio::osir::InfectionState::Count) &&
                        index < region * static_cast<size_t>(mio::osir::InfectionState::Count) +
                                    static_cast<size_t>(mio::osir::InfectionState::Count)) {
                        moment_results.get_last_value()[index] = 0;
                    }
                    // Moments
                    if (index >= m_moment_simulation.get_model().populations.get_num_compartments()) {
                        auto multiindex = m_moment_simulation.get_model().moments.unflatten_index(
                            index - m_moment_simulation.get_model().populations.get_num_compartments());
                        int order_region_to = std::accumulate(
                            multiindex.begin() + static_cast<size_t>(mio::osir::InfectionState::Count) * region,
                            multiindex.begin() + static_cast<size_t>(mio::osir::InfectionState::Count) * region +
                                static_cast<size_t>(mio::osir::InfectionState::Count),
                            0.0);
                        if (order_region_to >
                            0) { // Set moments that are fully or partially in stochastic region to zero
                            moment_results.get_last_value()[index] = 0;
                        }
                    }
                }
            }
        }

        // Update mean and moment ts of simulation set
        m_stochastic_simulation.get_mean().remove_last_time_point();
        m_stochastic_simulation.recalculate_last_mean();
        m_stochastic_simulation.get_moments().remove_last_time_point();
        m_stochastic_simulation.recalculate_last_moments();
    }

    void exchange_stochastic_to_moment()
    {
        auto moment_results            = m_moment_simulation.get_result();
        auto& stochastic_sims          = m_stochastic_simulation.get_simulations();
        auto& stochastic_results       = m_stochastic_simulation.get_result();
        auto& stochastic_model_means   = m_stochastic_simulation.get_mean();
        auto& stochastic_model_moments = m_stochastic_simulation.get_moments();
        for (size_t region_to : m_deterministic_regions) {
            bool region_exchanged = false;
            for (size_t state = 0; state < static_cast<size_t>(mio::osir::InfectionState::Count); ++state) {
                double mean =
                    stochastic_model_means
                        .get_last_value()[region_to * static_cast<size_t>(mio::osir::InfectionState::Count) + state];
                if (mean < 1.) { // No (or not enough) agents to exchange for this region and state
                    continue;
                }
                else {
                    region_exchanged = true;
                    // Add mean to moment model result
                    moment_results.get_last_value()[m_moment_simulation.get_model().populations.get_flat_index(
                        {mio::regions::Region(region_to), mio::osir::InfectionState(state)})] += mean;

                    // Reset results and simulations of simulation set
                    for (size_t sim = 0; sim < stochastic_sims.size(); ++sim) {
                        auto& sim_result = stochastic_sims[sim].get_result();
                        auto& sim_pop    = stochastic_sims[sim].get_model().populations;

                        // Reset last result in all Simulation objects
                        sim_result.get_last_value()[static_cast<size_t>(mio::osir::InfectionState::Count) * region_to +
                                                    state] = 0;
                        // Reset population in all Simulation objects
                        sim_pop[{mio::regions::Region(region_to), mio::osir::InfectionState(state)}] = 0;
                        // Reset last result in SimulationSet object
                        stochastic_sims[sim]
                            .get_last_value()[static_cast<size_t>(mio::osir::InfectionState::Count) * region_to +
                                              state] = 0;
                    }
                }
            }

            if (region_exchanged) {
                // Exchange moments
                auto moment_indices_stochastic_model = m_stochastic_simulation.get_moment_indices();
                for (size_t i = 0; i < moment_indices_stochastic_model.size(); ++i) {
                    int order = std::accumulate(moment_indices_stochastic_model[i].begin(),
                                                moment_indices_stochastic_model[i].end(), 0.0);
                    int order_region_to =
                        std::accumulate(moment_indices_stochastic_model[i].begin() +
                                            static_cast<size_t>(mio::osir::InfectionState::Count) * region_to,
                                        moment_indices_stochastic_model[i].begin() +
                                            static_cast<size_t>(mio::osir::InfectionState::Count) * region_to +
                                            static_cast<size_t>(mio::osir::InfectionState::Count),
                                        0.0);
                    if (order_region_to !=
                        0) { // If the moment has no indices in region_to, it does not have to be considered
                        if (order_region_to ==
                            order) { // Only moments that are fully in the region are considered (Variances)
                            auto index = m_moment_simulation.get_model().moments.flatten_index(
                                moment_indices_stochastic_model[i]);
                            // Sample new variance value
                            double var_stochastic = stochastic_model_moments.get_last_value()[i];
                            double var_moment =
                                moment_results.get_last_value()
                                    [m_moment_simulation.get_model().populations.get_num_compartments() + index];
                            double value =
                                std::max(0., mio::UniformDistribution<double>::get_instance()(
                                                 mio::thread_local_rng(),
                                                 std::pow(std::sqrt(var_moment) - std::sqrt(var_stochastic), 2.),
                                                 std::pow(std::sqrt(var_moment) + std::sqrt(var_stochastic), 2.)));
                            moment_results
                                .get_last_value()[m_moment_simulation.get_model().populations.get_num_compartments() +
                                                  index] = value;
                        }
                    }
                }
            }
        }

        // Reset and update mean and moment ts of stochastic model
        stochastic_model_means.remove_last_time_point();
        m_stochastic_simulation.recalculate_last_mean();
        stochastic_model_moments.remove_last_time_point();
        m_stochastic_simulation.recalculate_last_moments();
    }

    SMMSetSim initialize_stochastic_model(size_t num_runs)
    {
        assert(m_config->num_regions == num_regions);
        mio::smm::Model<ScalarType, num_regions, mio::osir::InfectionState> model =
            smm_helper::initialize_model<num_regions>(*m_config);
        return SMMSetSim(model, num_runs, m_config->t0, m_config->dt);
    }

    MomentSim initialize_deterministic_model(ClosureFunctionType closure_func, double min_step_size)
    {
        assert(m_config->num_regions == num_regions);
        mio::smm_moments::Model<num_regions, closure_order> model(closure_func);

        // Initialize model adoption and transition rates with zero as all regions start with the stochastic formulation
        // Set spatial transition rates
        for (auto& rate : m_config->transition_rates) {
            model.parameters.template get<mio::smm_moments::TransitionRate>()[{rate.status, rate.from, rate.to}] = 0.;
        }

        model.parameters.template get<mio::smm_moments::TransmissionRate>() = 0.0;
        model.parameters.template get<mio::smm_moments::RecoveryRate>()     = 0.0;

        MomentSim sim(model, m_config->t0, m_config->dt);
        if (min_step_size > 0) {
            sim.get_integrator_core().get_dt_min() = min_step_size;
        }
        return sim;
    }

    SMMSetSim m_stochastic_simulation; ///< SMM simulation containing stochastically modeled regions.
    MomentSim m_moment_simulation; ///< Moment simulation containing deterministically modeled regions.
    std::vector<size_t> m_stochastic_regions; ///< Regions which are currently modeled stochastically.
    std::vector<size_t> m_deterministic_regions; ///< Regions which are currently modeled deterministically.
    double m_t; ///< Current time step.
    double m_dt; ///< Step size.
    std::shared_ptr<Config::Config> m_config; ///< Config of the simulation.
};

} // namespace hybrid

} // namespace mio

#endif //MIO_SPATIAL_HYBRID_MODEL_H
