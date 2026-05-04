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
#include "memilio/data/analyze_result.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/random_number_generator.h"
#include "memilio/utils/time_series.h"
#include "ode_sir/infection_state.h"
#include "smm/parameters.h"
#include "smm/simulation_set.h"
#include "smm_moments/parameters.h"
#include "smm_moments/simulation.h"
#include "smm_moments/closure_functions.h"
#include "smm/model.h"
#include "simulations/hybrid_simulations/sir_metapop/config/config.h"
#include "simulations/hybrid_simulations/sir_metapop/library/smm_helper.h"
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <math.h>
#include <memory>
#include <numeric>
#include <utility>
#include <vector>
#include <unordered_set>

#include <boost/math/tools/roots.hpp>
#include <boost/math/distributions/normal.hpp>

namespace mio
{
namespace hybrid
{

template <size_t num_regions, size_t closure_order, typename ConfigType>
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
        const ConfigType& config, size_t num_runs, double min_step_size, double dt,
        ClosureFunctionType closure_func = &mio::smm_moments::truncation_closure<num_regions, closure_order>)
        : m_config(std::make_shared<ConfigType>(config))
        , m_stochastic_simulation(initialize_stochastic_model(num_runs))
        , m_moment_simulation(initialize_deterministic_model(closure_func, min_step_size))
        , m_stochastic_regions(num_regions)
        , m_deterministic_regions(0)
        , m_t(config.t0)
        , m_dt(dt)
        , m_model_used(Eigen::Index(num_regions))
    {
        // Initially all regions are modeled stochastically
        std::iota(m_stochastic_regions.begin(), m_stochastic_regions.end(), 0);
    }

    /**
     * @brief Advance spatial-hybrid simulation.
     * @param[in] tmax Simulation end time point.
     * @param[in] condition Switching condition for determining when to switch a region between stochastic and deterministic modeling.
     */
    template <typename Condition>
    void advance(double tmax, Condition&& condition, bool log_model_used = false)
    {
        while (m_t < tmax) {
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
            { // Log information which model is used at current time point
                if (log_model_used) {
                    Eigen::VectorXi zeros = Eigen::VectorXi::Constant(num_regions, 0);
                    m_model_used.add_time_point(m_t, zeros);
                    for (size_t det_region : m_deterministic_regions) {
                        m_model_used.get_last_value()[det_region] = 1;
                    }
                }
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

    /**
     * @brief Get mean time series of stochastic simulation.
     * The time series is interpolated to the given interpolation step size.
     * @param[in] interpolation_step_size Step size for interpolation of mean time series.
     */
    TimeSeries<double> get_stochastic_mean(double interpolation_step_size)
    {
        auto mean_ts  = m_stochastic_simulation.get_mean();
        int num_steps = static_cast<int>((mean_ts.get_last_time() - mean_ts.get_time(0)) / interpolation_step_size) + 1;
        std::vector<double> interpolation_tps(num_steps);

        for (int i = 0; i < num_steps; ++i) {
            interpolation_tps[i] = i * interpolation_step_size;
        }
        return interpolate_simulation_result(mean_ts, interpolation_tps);
    }

    /**
     * @brief Get moment time series of stochastic simulation.
     * The time series is interpolated to the given interpolation step size.
     * @param[in] interpolation_step_size Step size for interpolation of moment time series.
     */
    std::pair<TimeSeries<double>, std::vector<std::string>> get_stochastic_moments(double interpolation_step_size)
    {
        auto moment_ts = m_stochastic_simulation.get_moments();
        int num_steps =
            static_cast<int>((moment_ts.get_last_time() - moment_ts.get_time(0)) / interpolation_step_size) + 1;
        std::vector<double> interpolation_tps(num_steps);

        for (int i = 0; i < num_steps; ++i) {
            interpolation_tps[i] = i * interpolation_step_size;
        }
        return std::make_pair(interpolate_simulation_result(moment_ts, interpolation_tps),
                              m_stochastic_simulation.get_moment_names());
    }

    /**
     * @brief Get mean time series of moment simulation.
     * The time series is interpolated to the given interpolation step size.
     * @param[in] interpolation_step_size Step size for interpolation of mean time series.
     */
    TimeSeries<double> get_deterministic_mean(double interpolation_step_size)
    {
        auto mean_ts  = m_moment_simulation.get_expected_values_time_series();
        int num_steps = static_cast<int>((mean_ts.get_last_time() - mean_ts.get_time(0)) / interpolation_step_size) + 1;
        std::vector<double> interpolation_tps(num_steps);

        for (int i = 0; i < num_steps; ++i) {
            interpolation_tps[i] = i * interpolation_step_size;
        }
        return interpolate_simulation_result(mean_ts, interpolation_tps);
    }

    /**
     * @brief Get moment time series of moment simulation.
     * The time series is interpolated to the given interpolation step size.
     * @param[in] interpolation_step_size Step size for interpolation of moment time series.
     */
    std::pair<TimeSeries<double>, std::vector<std::string>> get_deterministic_moments(double interpolation_step_size)
    {
        auto moment_ts = m_moment_simulation.get_moment_time_series(closure_order);
        int num_steps  = static_cast<int>((moment_ts.first.get_last_time() - moment_ts.first.get_time(0)) /
                                         interpolation_step_size) +
                        1;
        std::vector<double> interpolation_tps(num_steps);

        for (int i = 0; i < num_steps; ++i) {
            interpolation_tps[i] = i * interpolation_step_size;
        }
        return std::make_pair(interpolate_simulation_result(moment_ts.first, interpolation_tps), moment_ts.second);
    }

    TimeSeries<double> get_joint_mean(double interpolation_step_size)
    {
        return merge_time_series(get_stochastic_mean(interpolation_step_size),
                                 get_deterministic_mean(interpolation_step_size), true)
            .value();
    }

    std::pair<TimeSeries<double>, std::vector<std::string>> get_joint_moments(double interpolation_step_size)
    {
        auto stochastic_moments    = get_stochastic_moments(interpolation_step_size);
        auto deterministic_moments = get_deterministic_moments(interpolation_step_size);
        assert(stochastic_moments.second.size() == deterministic_moments.second.size());
        for (size_t i = 0; i < stochastic_moments.second.size(); ++i) {
            assert(stochastic_moments.second[i] == deterministic_moments.second[i]);
        }
        return std::make_pair(merge_time_series(stochastic_moments.first, deterministic_moments.first, true).value(),
                              stochastic_moments.second);
    }

    TimeSeries<int> get_model_used_ts()
    {
        return m_model_used;
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
                if (std::accumulate(moment_index.begin(), moment_index.end(), 0.0) > 0) {
                    assert(moment_model_result.get_last_value()[flat_index] ==
                           0.); // Check that region is not already modeled in moment model
                    moment_model_result.get_last_value()[flat_index] = stochastic_model_moments.get_last_value()[i];
                }
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

#ifdef MEMILIO_ENABLE_OPENMP
#pragma omp parallel for
#endif
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

        // Set means and moments in corresponding regions to zero
        for (size_t region_to : regions_to_switch) {
            m_stochastic_simulation.remove_means_for_region(region_to);
            m_stochastic_simulation.remove_moments_for_region(region_to);
        }
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
            if (typeid(*m_config) == typeid(Config::sirs::Config)) {
                // Set immunity loss rate for sirs
                m_moment_simulation.get_model()
                    .parameters.template get<mio::smm_moments::ImmunityLossRate>()[mio::regions::Region(region)] =
                    m_config->nu;
            }

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
#ifdef MEMILIO_ENABLE_OPENMP
#pragma omp parallel for
#endif
            for (size_t sim = 0; sim < m_stochastic_simulation.get_simulations().size(); ++sim) {
                auto& sim_model = m_stochastic_simulation.get_simulations()[sim].get_model();
                for (auto& rate :
                     sim_model.parameters.template get<mio::smm::AdoptionRates<double, mio::osir::InfectionState>>()) {
                    if (rate.region == mio::regions::Region(region)) {
                        rate.factor = 0.;
                    }
                }

                for (auto& rate : sim_model.parameters
                                      .template get<mio::smm::TransitionRates<double, mio::osir::InfectionState>>()) {
                    if (rate.from == mio::regions::Region(region)) {
                        rate.factor = 0.;
                    }
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

                // Get current mean and variance from moment model for region and state
                double mean = moment_results[region * static_cast<size_t>(mio::osir::InfectionState::Count) + state];
                std::array<int, static_cast<size_t>(mio::osir::InfectionState::Count) * num_regions> indices_var;
                indices_var.fill(0);
                indices_var[region * static_cast<size_t>(mio::osir::InfectionState::Count) + state] = 2;
                double std =
                    std::sqrt(moment_results[m_moment_simulation.get_model().populations.get_num_compartments() +
                                             m_moment_simulation.get_model().moments.flatten_index(indices_var)]);
                if (mean - 3 * std < 0) {
                    auto new_params = calculate_trunc_normal_params(0., mean, std);
                    double l_old    = lambda(-mean / std);
                    double old_m    = mean + std * l_old;
                    double old_s    = std * std * (1 + (-mean / std) * l_old - l_old * l_old);
                    double l_new    = lambda(-new_params.first / new_params.second);
                    double new_m    = new_params.first + new_params.second * l_new;
                    double new_s    = new_params.second * new_params.second *
                                   (1 + (-new_params.first / new_params.second) * l_new - l_new * l_new);
                    if (std::abs(mean - old_m) > std::abs(mean - new_m)) {
                        mean = new_params.first;
                    }
                    if (std::abs(std - old_s) > std::abs(std - new_s)) {
                        std = new_params.second;
                    }
                }
                double rounding_err = 0.;
#ifdef MEMILIO_ENABLE_OPENMP
#pragma omp parallel for
#endif
                for (size_t sim = 0; sim < stochastic_sims.size(); ++sim) {
                    auto& sim_result = stochastic_sims[sim].get_result();
                    auto& sim_pop    = stochastic_sims[sim].get_model().populations;

                    // Sample number of agents for stochastic trajectory normally distributed
                    double sim_value =
                        sample_truncated_normal(stochastic_sims[sim].get_model().get_rng(), mean, std, 0.);
                    while (sim_value < 0) {
                        sim_value = sample_truncated_normal(stochastic_sims[sim].get_model().get_rng(), mean, std, 0.);
                        // sim_value = mio::NormalDistribution<double>::get_instance()(
                        //     stochastic_sims[sim].get_model().get_rng(), mean, std);
                    }
                    if (correct_for_rounding_error && rounding_err != 0.) {
                        sim_value += rounding_err;
                        rounding_err = 0.;
                    }
                    double sim_value_rounded = std::max(0., std::round(sim_value));
                    rounding_err += sim_value - sim_value_rounded;
                    // Set population in Simulation object
                    sim_pop[{mio::regions::Region(region), mio::osir::InfectionState(state)}] = sim_value_rounded;

                    // Set sampled value in SMMSet result object
                    assert(m_moment_simulation.get_result().get_last_time() == stochastic_results[sim].get_last_time());
                    stochastic_results[sim]
                        .get_last_value()[region * static_cast<size_t>(mio::osir::InfectionState::Count) + state] =
                        sim_value_rounded;
                    // Set sampled value in Simulation object
                    assert(m_moment_simulation.get_result().get_last_time() == sim_result.get_last_time());
                    sim_result
                        .get_last_value()[region * static_cast<size_t>(mio::osir::InfectionState::Count) + state] =
                        sim_value_rounded;
                }
            }
        }

        //Reset and update mean and moment ts of simulation set
        auto& stochastic_model_means   = m_stochastic_simulation.get_mean();
        auto& stochastic_model_moments = m_stochastic_simulation.get_moments();
        stochastic_model_means.remove_last_time_point();
        stochastic_model_moments.remove_last_time_point();
        m_stochastic_simulation.recalculate_last_means_and_moments();
    }

    void set_rates_in_stochastic_model(const std::unordered_set<size_t>& regions_to_switch)
    {
        if (regions_to_switch.empty()) {
            return;
        }
        for (size_t region : regions_to_switch) {
#ifdef MEMILIO_ENABLE_OPENMP
#pragma omp parallel for
#endif
            for (size_t sim = 0; sim < m_stochastic_simulation.get_simulations().size(); ++sim) {
                auto& sim_model = m_stochastic_simulation.get_simulations()[sim].get_model();
                for (auto& rate :
                     sim_model.parameters.template get<mio::smm::AdoptionRates<double, mio::osir::InfectionState>>()) {
                    if (rate.region == mio::regions::Region(region)) {
                        if (rate.from == mio::osir::InfectionState::Susceptible) {
                            rate.factor = m_config->lambdas[region];
                        }
                        if (rate.from == mio::osir::InfectionState::Infected) {
                            rate.factor = m_config->gamma;
                        }
                        if (typeid(*m_config) == typeid(Config::sirs::Config)) {
                            if (rate.from == mio::osir::InfectionState::Recovered) {
                                rate.factor = m_config->nu;
                            }
                        }
                    }
                }

                for (auto& rate : sim_model.parameters
                                      .template get<mio::smm::TransitionRates<double, mio::osir::InfectionState>>()) {
                    if (rate.from == mio::regions::Region(region)) {
                        auto config_rate_it = std::find_if(
                            m_config->transition_rates.begin(), m_config->transition_rates.end(), [&](const auto& r) {
                                return r.status == rate.status && r.from == rate.from && r.to == rate.to;
                            });
                        rate.factor = config_rate_it->factor;
                    }
                }
            }
        }
    }

    void remove_regions_from_moment_model(const std::unordered_set<size_t>& regions_to_switch)
    {
        if (regions_to_switch.empty()) {
            return;
        }
        auto& moment_result = m_moment_simulation.get_result();
        for (size_t region : regions_to_switch) {
            for (size_t index = 0; index < static_cast<size_t>(moment_result.get_num_elements()); ++index) {
                // Mean entries for corresponding regions are set to zero
                if (index >= region * static_cast<size_t>(mio::osir::InfectionState::Count) &&
                    index < region * static_cast<size_t>(mio::osir::InfectionState::Count) +
                                static_cast<size_t>(mio::osir::InfectionState::Count)) {
                    moment_result.get_last_value()[index]                                                          = 0;
                    m_moment_simulation.get_model().populations[{
                        mio::regions::Region(region),
                        mio::osir::InfectionState(index % static_cast<size_t>(mio::osir::InfectionState::Count))}] = 0;
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
            if (typeid(*m_config) == typeid(Config::sirs::Config)) {
                // Set immunity loss rate for sirs
                m_moment_simulation.get_model()
                    .parameters.template get<mio::smm_moments::ImmunityLossRate>()[mio::regions::Region(region)] = 0.;
            }

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
        auto& moment_results     = m_moment_simulation.get_result();
        auto& stochastic_sims    = m_stochastic_simulation.get_simulations();
        auto& stochastic_results = m_stochastic_simulation.get_result();
        bool recalc_stoch_means  = false;
        for (size_t region : m_stochastic_regions) {
            bool region_exchanged = false;
            for (size_t state = 0; state < static_cast<size_t>(mio::osir::InfectionState::Count); ++state) {
                double mean =
                    moment_results
                        .get_last_value()[region * static_cast<size_t>(mio::osir::InfectionState::Count) + state];
                if (mean <= 0.) { // No agents to exchange for this region and state
                    continue;
                }
                else {
                    recalc_stoch_means = true;
                    region_exchanged   = true;
                    std::array<int, static_cast<size_t>(mio::osir::InfectionState::Count) * num_regions> indices_var;
                    indices_var.fill(0);
                    indices_var[region * static_cast<size_t>(mio::osir::InfectionState::Count) + state] = 2;
                    double std                                                                          = std::sqrt(
                        moment_results
                            .get_last_value()[m_moment_simulation.get_model().populations.get_num_compartments() +
                                              m_moment_simulation.get_model().moments.flatten_index(indices_var)]);
                    if (mean - 3 * std < 0) {
                        auto new_params = calculate_trunc_normal_params(0., mean, std);
                        double l_old    = lambda(-mean / std);
                        double old_m    = mean + std * l_old;
                        double old_s    = std * std * (1 + (-mean / std) * l_old - l_old * l_old);
                        double l_new    = lambda(-new_params.first / new_params.second);
                        double new_m    = new_params.first + new_params.second * l_new;
                        double new_s    = new_params.second * new_params.second *
                                       (1 + (-new_params.first / new_params.second) * l_new - l_new * l_new);
                        if (std::abs(mean - old_m) > std::abs(mean - new_m)) {
                            mean = new_params.first;
                        }
                        if (std::abs(std - old_s) > std::abs(std - new_s)) {
                            std = new_params.second;
                        }
                    }
                    double rounding_err = 0.;
#ifdef MEMILIO_ENABLE_OPENMP
#pragma omp parallel for
#endif
                    for (size_t sim = 0; sim < stochastic_sims.size(); ++sim) {
                        auto& sim_result = stochastic_sims[sim].get_result();
                        auto& sim_pop    = stochastic_sims[sim].get_model().populations;

                        // Sample number of incoming agents for simulation
                        double sim_value =
                            sample_truncated_normal(stochastic_sims[sim].get_model().get_rng(), mean, std, 0.);
                        while (sim_value < 0) {
                            sim_value =
                                sample_truncated_normal(stochastic_sims[sim].get_model().get_rng(), mean, std, 0.);
                        }
                        if (correct_for_rounding_error && rounding_err != 0.) {
                            sim_value += rounding_err;
                            rounding_err = 0.;
                        }
                        double sim_value_rounded = std::max(0., std::round(sim_value));
                        rounding_err += sim_value - sim_value_rounded;

                        // Assert time last time point of moment and stochastic simulation is the same
                        assert(m_moment_simulation.get_result().get_last_time() == sim_result.get_last_time());
                        assert(m_moment_simulation.get_result().get_last_time() ==
                               stochastic_results[sim].get_last_time());
                        // Add agents to last result in all Simulation objects
                        sim_result
                            .get_last_value()[region * static_cast<size_t>(mio::osir::InfectionState::Count) + state] +=
                            sim_value_rounded;
                        // Add agents to population in all Simulation objects
                        sim_pop[{mio::regions::Region(region), mio::osir::InfectionState(state)}] += sim_value_rounded;
                        // Add agents to last result in SimulationSet object
                        stochastic_results[sim]
                            .get_last_value()[static_cast<size_t>(mio::osir::InfectionState::Count) * region + state] +=
                            sim_value_rounded;
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

        if (recalc_stoch_means) { // Update mean and moment ts of simulation set
            m_stochastic_simulation.get_mean().remove_last_time_point();
            m_stochastic_simulation.get_moments().remove_last_time_point();
            m_stochastic_simulation.recalculate_last_means_and_moments();
        }
    }

    void exchange_stochastic_to_moment()
    {
        auto& moment_results                  = m_moment_simulation.get_result();
        auto& stochastic_sims                 = m_stochastic_simulation.get_simulations();
        auto& stochastic_model_means          = m_stochastic_simulation.get_mean();
        auto& stochastic_model_moments        = m_stochastic_simulation.get_moments();
        auto& moment_indices_stochastic_model = m_stochastic_simulation.get_moment_indices();
        auto& moment_orders_stochastic_model  = m_stochastic_simulation.get_moment_orders();
        bool recalc_stoch_mean                = false;
        for (size_t region_to : m_deterministic_regions) {
            bool region_exchanged = false;
            for (size_t state = 0; state < static_cast<size_t>(mio::osir::InfectionState::Count); ++state) {
                double mean =
                    stochastic_model_means
                        .get_last_value()[region_to * static_cast<size_t>(mio::osir::InfectionState::Count) + state];
                if (mean <= 0.) { // No agents to exchange for this region and state
                    continue;
                }
                else {
                    recalc_stoch_mean = true;
                    region_exchanged  = true;
                    // Add mean to moment model result
                    moment_results.get_last_value()[m_moment_simulation.get_model().populations.get_flat_index(
                        {mio::regions::Region(region_to), mio::osir::InfectionState(state)})] += mean;

#ifdef MEMILIO_ENABLE_OPENMP
#pragma omp parallel for
#endif
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
                        sim_result.get_last_value()[static_cast<size_t>(mio::osir::InfectionState::Count) * region_to +
                                                    state] = 0;
                    }
                }
            }

            if (region_exchanged) {
                // Exchange moments
                for (size_t i = 0; i < moment_indices_stochastic_model.size(); ++i) {
                    int order = moment_orders_stochastic_model[i];
                    int order_region_to =
                        std::accumulate(moment_indices_stochastic_model[i].begin() +
                                            static_cast<size_t>(mio::osir::InfectionState::Count) * region_to,
                                        moment_indices_stochastic_model[i].begin() +
                                            static_cast<size_t>(mio::osir::InfectionState::Count) * region_to +
                                            static_cast<size_t>(mio::osir::InfectionState::Count),
                                        0.0);
                    if (order_region_to !=
                        0) { // If the moment has no indices in region_to, it does not have to be considered
                        if (order_region_to == order &&
                            (std::count(moment_indices_stochastic_model[i].begin(),
                                        moment_indices_stochastic_model[i].end(), 2) == 1 &&
                             std::count(moment_indices_stochastic_model[i].begin(),
                                        moment_indices_stochastic_model[i].end(),
                                        0) == num_regions * static_cast<size_t>(mio::osir::InfectionState::Count) -
                                                  1)) { // Only variances that are fully in the region are considered
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

        if (recalc_stoch_mean) {
            for (size_t region_to : m_deterministic_regions) {
                // Remove means and moments in deterministic regions
                m_stochastic_simulation.remove_means_for_region(region_to);
                m_stochastic_simulation.remove_moments_for_region(region_to);
            }
        }
    }

    SMMSetSim initialize_stochastic_model(size_t num_runs)
    {
        assert(m_config->num_regions == num_regions);
        mio::smm::Model<ScalarType, num_regions, mio::osir::InfectionState> model =
            smm_helper::initialize_model<num_regions>(*m_config);
        return SMMSetSim(num_runs, model, m_config->t0, m_config->dt);
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
        model.parameters.template get<mio::smm_moments::ImmunityLossRate>() = 0.0;

        MomentSim sim(model, m_config->t0, m_config->dt);
        if (min_step_size > 0) {
            sim.get_integrator_core().get_dt_min() = min_step_size;
        }
        sim.get_integrator_core().get_dt_max() = m_config->dt;
        return sim;
    }

    double Phi(double x)
    {
        return 0.5 * (1. + std::erf(x / std::sqrt(2.)));
    }

    double phi(double x)
    {
        static const double INV_SQRT_2PI = 1.0 / std::sqrt(2.0 * M_PI);
        return INV_SQRT_2PI * std::exp(-0.5 * x * x);
    }

    double lambda(double alpha)
    {
        return phi(alpha) / (1. - Phi(alpha));
    }

    double V(double alpha)
    {
        double l = lambda(alpha);
        return 1 + alpha * l - l * l;
    }

    double f(double alpha, double z)
    {
        double l = lambda(alpha);
        double v = V(alpha);

        if (v < 0) {
            v = 0;
        }
        return (l - alpha) / std::sqrt(V(alpha)) - z;
    }

    double solve_for_alpha(double m, double a, double s)
    {
        auto f = [m, a, s, this](double alpha) {
            double l = lambda(alpha);
            double v = V(alpha);
            if (v <= 0) {
                v = 1;
            }
            return (l - alpha) / std::sqrt(v) - (m - a) / s;
        };
        double lo = -1.0;
        double hi = 1.0;

        while (f(lo) * f(hi) > 0) {
            lo *= 2.;
            hi *= 2.;
        }

        auto tol = [](double x1, double x2) {
            return std::abs(x1 - x2) < 1e-8;
        };

        std::uintmax_t max_iter = 100;

        std::pair<double, double> result = boost::math::tools::toms748_solve(f, lo, hi, tol, max_iter);
        return 0.5 * (result.first + result.second);
    }

    std::pair<double, double> calculate_trunc_normal_params(double a, double m, double s)
    {
        double alpha = solve_for_alpha(m, a, s);

        double sigma0 = s / std::sqrt(V(alpha));
        double mu0    = a - alpha * sigma0;
        return std::make_pair(mu0, sigma0);
    }

    double sample_truncated_normal(mio::RandomNumberGenerator& rng, double mu0, double sigma0, double a)
    {
        boost::math::normal dist;

        double alpha = (a - mu0) / sigma0;

        double Phi_alpha = cdf(dist, alpha);

        double u = mio::UniformDistribution<ScalarType>::get_instance()(rng, 0.0, 1.0);

        double u_prime = Phi_alpha + u * (1.0 - Phi_alpha);

        const double eps = 1e-15;

        if (u_prime >= 1.0)
            u_prime = 1.0 - eps;
        if (u_prime <= 0.0)
            u_prime = eps;

        double z = quantile(dist, u_prime);

        return mu0 + sigma0 * z;
    }

    std::shared_ptr<ConfigType> m_config; ///< Config of the simulation.
    SMMSetSim m_stochastic_simulation; ///< SMM simulation containing stochastically modeled regions.
    MomentSim m_moment_simulation; ///< Moment simulation containing deterministically modeled regions.
    std::vector<size_t> m_stochastic_regions; ///< Regions which are currently modeled stochastically.
    std::vector<size_t> m_deterministic_regions; ///< Regions which are currently modeled deterministically.
    double m_t; ///< Current time step.
    double m_dt; ///< Step size.
    TimeSeries<int>
        m_model_used; ///< Time series which indicates which model is used for each region at each time step (0: stochastic, 1: moment).

    const bool correct_for_rounding_error = true;
};

} // namespace hybrid

} // namespace mio

#endif //MIO_SPATIAL_HYBRID_MODEL_H
