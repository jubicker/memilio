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
#ifndef CONDITION_FUNCTIONS_H
#define CONDITION_FUNCTIONS_H

#include "simulations/hybrid_simulations/sir_metapop/config/config.h"
#include "smm/simulation_set.h"
#include "smm_moments/simulation.h"
#include <cstddef>
#include <memory>
#include <vector>

namespace mio
{

namespace hybrid
{

template <size_t num_regions, size_t closure_order>
class SwitchingCondition
{
    inline static double m_rel_switch_threshold            = 0.1;
    inline static double m_absolute_switch_threshold       = 100.;
    inline static double m_mean_stddev_relation            = 0.6;
    inline static double m_var_gradient_threshold          = 0;
    inline static double m_mean_gradient_threshold          = 0;
    inline static double m_timepoint_threshold             = 0;
    inline static std::unique_ptr<Config::Config> m_config = nullptr;

public:
    SwitchingCondition() = default;

    static void set_rel_switch_threshold(double value)
    {
        m_rel_switch_threshold = value;
    }

    static void set_config(const Config::Config& config)
    {
        m_config = std::make_unique<Config::Config>(config);
    }

    static void set_absolute_switch_threshold(double value)
    {
        m_absolute_switch_threshold = value;
    }

    static void set_mean_stddev_relation(double value)
    {
        m_mean_stddev_relation = value;
    }

    static void set_var_gradient_threshold(double value)
    {
        m_var_gradient_threshold = value;
    }

    static void set_mean_gradient_threshold(double value)
    {
        m_mean_gradient_threshold = value;
    }

    static void set_timepoint_threshold(double value)
    {
        m_timepoint_threshold = value;
    }

    static bool rel_threshold_condition(std::vector<double>& result_smm, std::vector<double>& result_moments,
                                        bool smm_used)
    {
        double total_population =
            std::accumulate(m_config->total_populations.begin(), m_config->total_populations.end(), 0.);
        if (smm_used) {
            double total_infected = 0;
            for (size_t r = 0; r < m_config->num_regions; ++r) {
                total_infected +=
                    result_smm[r * (int)mio::osir::InfectionState::Count + (int)mio::osir::InfectionState::Infected];
            }
            if ((total_infected > m_rel_switch_threshold * total_population) || (total_infected < 1)) {
                return true;
            }
        }
        else { // moment model currently used
            double total_infected = 0;
            for (size_t r = 0; r < m_config->num_regions; ++r) {
                total_infected += result_moments[r * (int)mio::osir::InfectionState::Count +
                                                 (int)mio::osir::InfectionState::Infected];
            }
            if ((total_infected < m_rel_switch_threshold * total_population) && (total_infected >= 1)) {
                return true;
            }
        }
        return false;
    }

    static bool abs_threshold_condition(std::vector<double>& result_smm, std::vector<double>& result_moments,
                                        bool smm_used)
    {
        if (smm_used) {
            double total_infected = 0;
            for (size_t r = 0; r < m_config->num_regions; ++r) {
                total_infected +=
                    result_smm[r * (int)mio::osir::InfectionState::Count + (int)mio::osir::InfectionState::Infected];
            }
            if ((total_infected > m_absolute_switch_threshold) || (total_infected < 1)) {
                return true;
            }
        }
        else { // moment model currently used
            double total_infected = 0;
            for (size_t r = 0; r < m_config->num_regions; ++r) {
                total_infected += result_moments[r * (int)mio::osir::InfectionState::Count +
                                                 (int)mio::osir::InfectionState::Infected];
            }
            if ((total_infected < m_absolute_switch_threshold) && (total_infected >= 1)) {
                return true;
            }
        }
        return false;
    }

    static bool var_gradient_condition(std::vector<double>& result_smm, std::vector<double>& result_moments,
                                       bool smm_used)
    {
        if (smm_used) {
            for (size_t r = 0; r < m_config->num_regions; ++r) {
                auto var_infected_gradient =
                    result_smm[r * (int)mio::osir::InfectionState::Count + (int)mio::osir::InfectionState::Infected];

                if (var_infected_gradient < m_var_gradient_threshold) {
                    return true;
                }
            }
        }
        else { // moment model currently used
            for (size_t r = 0; r < m_config->num_regions; ++r) {
                auto var_infected_gradient = result_moments[r * (int)mio::osir::InfectionState::Count +
                                                            (int)mio::osir::InfectionState::Infected];

                if (var_infected_gradient < m_var_gradient_threshold) {
                    return false;
                }
            }
            return true;
        }
        return false;
    }

    static bool mean_stddev_relation_condition(std::vector<double>& result_smm, std::vector<double>& result_moments,
                                               bool smm_used)
    {
        if (smm_used) {
            for (size_t r = 0; r < m_config->num_regions; ++r) {
                auto relation_infected =
                    result_smm[r * (int)mio::osir::InfectionState::Count + (int)mio::osir::InfectionState::Infected];
                if (relation_infected > m_mean_stddev_relation) {
                    return false;
                }
                return true;
            }
        }
        else { // moment model currently used
            for (size_t r = 0; r < m_config->num_regions; ++r) {
                auto relation_infected = result_moments[r * (int)mio::osir::InfectionState::Count +
                                                        (int)mio::osir::InfectionState::Infected];
                if (relation_infected > m_mean_stddev_relation) {
                    return true;
                }
            }
        }
        return false;
    }

    static bool region_index_condition(
        mio::smm::SimulationSet<num_regions, mio::osir::InfectionState, closure_order>& /*stochastic_model*/,
        mio::smm_moments::Simulation<num_regions, closure_order>& /*deterministic_model*/, bool stochastic_used,
        size_t region)
    {
        // Regions with even index are modeled deterministically, regions with odd index are modeled stochastically
        return stochastic_used ? (region % 2 == 0) : (region % 2 == 1);
    }

    static bool reversed_region_index_condition(
        mio::smm::SimulationSet<num_regions, mio::osir::InfectionState, closure_order>& /*stochastic_model*/,
        mio::smm_moments::Simulation<num_regions, closure_order>& /*deterministic_model*/, bool stochastic_used,
        size_t region)
    {
        // Regions with odd index are modeled deterministically, regions with even index are modeled stochastically
        return stochastic_used ? (region % 2 == 1) : (region % 2 == 0);
    }

    static bool mean_stddev_relation_region(
        mio::smm::SimulationSet<num_regions, mio::osir::InfectionState, closure_order>& stochastic_model,
        mio::smm_moments::Simulation<num_regions, closure_order>& deterministic_model, bool stochastic_used,
        size_t region)
    {
        if (stochastic_used) {
            auto relations = current_smm_relations(stochastic_model, 0.);
            if ((relations[region * (int)mio::osir::InfectionState::Count + (int)mio::osir::InfectionState::Infected] >
                0) && (relations[region * (int)mio::osir::InfectionState::Count + (int)mio::osir::InfectionState::Infected] <
                m_mean_stddev_relation)) {
                return true;
            }
            return false;
        }
        else {
            auto relations = current_moment_relations(deterministic_model, 0.);
            if ((relations[region * (int)mio::osir::InfectionState::Count + (int)mio::osir::InfectionState::Infected] <=
                0) || (relations[region * (int)mio::osir::InfectionState::Count + (int)mio::osir::InfectionState::Infected] >
                m_mean_stddev_relation)) {
                return true;
            }
            return false;
        }
    }

    static bool combined_relation_var_gradient_condition_region(
        mio::smm::SimulationSet<num_regions, mio::osir::InfectionState, closure_order>& stochastic_model,
        mio::smm_moments::Simulation<num_regions, closure_order>& deterministic_model, bool stochastic_used,
        size_t region)
    {
        if (stochastic_used) {
            auto relations     = current_smm_relations(stochastic_model, 0.);
            auto var_gradients = current_smm_var_gradients(stochastic_model, 0.);
            if ((var_gradients[region * (int)mio::osir::InfectionState::Count +
                               (int)mio::osir::InfectionState::Infected] <= m_var_gradient_threshold) &&
                (relations[region * (int)mio::osir::InfectionState::Count + (int)mio::osir::InfectionState::Infected] <
                 m_mean_stddev_relation)) {
                return true;
            }
            return false;
        }
        else {
            auto relations     = current_moment_relations(deterministic_model, 0.);
            auto var_gradients = current_moment_var_gradients(deterministic_model, 0.);
            if ((var_gradients[region * (int)mio::osir::InfectionState::Count +
                               (int)mio::osir::InfectionState::Infected] > m_var_gradient_threshold + 10) &&
                relations[region * (int)mio::osir::InfectionState::Count + (int)mio::osir::InfectionState::Infected] >
                    1.1 * m_mean_stddev_relation) {
                return true;
            }
            return false;
        }
    }

    static bool
    fixed_tp_region(mio::smm::SimulationSet<num_regions, mio::osir::InfectionState, closure_order>& stochastic_model,
                    mio::smm_moments::Simulation<num_regions, closure_order>& /*deterministic_model*/,
                    bool stochastic_used, size_t /*region*/)
    {
        if (stochastic_used) {
            if (stochastic_model.get_mean().get_last_time() > m_timepoint_threshold) {
                return true;
            }
        }
        return false;
    }

    static bool abs_threshold_region(
        mio::smm::SimulationSet<num_regions, mio::osir::InfectionState, closure_order>& stochastic_model,
        mio::smm_moments::Simulation<num_regions, closure_order>& deterministic_model, bool stochastic_used,
        size_t region)
    {
        if (stochastic_used) {
            auto current_means = current_smm_means(stochastic_model, 0.);
            if ((current_means[region * (int)mio::osir::InfectionState::Count +
                               (int)mio::osir::InfectionState::Infected] > m_absolute_switch_threshold) ||
                (current_means[region * (int)mio::osir::InfectionState::Count +
                               (int)mio::osir::InfectionState::Infected] < 1)) {
                return true;
            }
        }
        else { // moment model currently used
            auto current_means = current_moment_means(deterministic_model, 0.);
            if ((current_means[region * (int)mio::osir::InfectionState::Count +
                               (int)mio::osir::InfectionState::Infected] < m_absolute_switch_threshold) &&
                (current_means[region * (int)mio::osir::InfectionState::Count +
                               (int)mio::osir::InfectionState::Infected] >= 1)) {
                return true;
            }
        }
        return false;
    }

    static bool R0_relation_condition_region(
        mio::smm::SimulationSet<num_regions, mio::osir::InfectionState, closure_order>& stochastic_model,
        mio::smm_moments::Simulation<num_regions, closure_order>& deterministic_model, bool stochastic_used,
        size_t region)
    {
        if (stochastic_used) {
            auto relations     = current_smm_relations(stochastic_model, 0.);
            auto current_means = current_smm_means(stochastic_model, 0.);
            double R0          = m_config->lambdas[region] / m_config->gamma *
                        current_means[region * (int)mio::osir::InfectionState::Count +
                                      (int)mio::osir::InfectionState::Susceptible];
            if ((R0 < 1.0 || R0 > 3.0) &&
                (relations[region * (int)mio::osir::InfectionState::Count + (int)mio::osir::InfectionState::Infected] <
                 m_mean_stddev_relation)) {
                return true;
            }
            return false;
        }
        else {
            auto relations     = current_moment_relations(deterministic_model, 0.);
            auto current_means = current_moment_means(deterministic_model, 0.);
            double R0          = m_config->lambdas[region] / m_config->gamma *
                        current_means[region * (int)mio::osir::InfectionState::Count +
                                      (int)mio::osir::InfectionState::Susceptible];
            if ((R0 >= 1.0 && R0 <= 3.0) &&
                relations[region * (int)mio::osir::InfectionState::Count + (int)mio::osir::InfectionState::Infected] >
                    1.1 * m_mean_stddev_relation) {
                return true;
            }
            return false;
        }
    }

    static bool R0_region(
        mio::smm::SimulationSet<num_regions, mio::osir::InfectionState, closure_order>& stochastic_model,
        mio::smm_moments::Simulation<num_regions, closure_order>& deterministic_model, bool stochastic_used,
        size_t region)
    {
        if (stochastic_used) {
            auto current_means = current_smm_means(stochastic_model, 0.);
            double R0          = m_config->lambdas[region] / m_config->gamma *
                        current_means[region * (int)mio::osir::InfectionState::Count +
                                      (int)mio::osir::InfectionState::Susceptible];
            if (R0 < 1.0) {
                return true;
            }
            return false;
        }
        else {
            auto current_means = current_moment_means(deterministic_model, 0.);
            double R0          = m_config->lambdas[region] / m_config->gamma *
                        current_means[region * (int)mio::osir::InfectionState::Count +
                                      (int)mio::osir::InfectionState::Susceptible];
            if (R0 > 1.0) {
                return true;
            }
            return false;
        }
    }
    static bool mean_gradient_region(
        mio::smm::SimulationSet<num_regions, mio::osir::InfectionState, closure_order>& stochastic_model,
        mio::smm_moments::Simulation<num_regions, closure_order>& deterministic_model, bool stochastic_used,
        size_t region)
    {
        if (stochastic_used) {
            auto var_gradients = current_smm_mean_gradients(stochastic_model, 0.);
            if ((var_gradients[region * (int)mio::osir::InfectionState::Count +
                               (int)mio::osir::InfectionState::Infected] <= m_mean_gradient_threshold)) {
                return true;
            }
            return false;
        }
        else {
            auto var_gradients = current_moment_mean_gradients(deterministic_model, 0.);
            if ((var_gradients[region * (int)mio::osir::InfectionState::Count +
                               (int)mio::osir::InfectionState::Infected] > m_mean_gradient_threshold)) {
                return true;
            }
            return false;
        }
    }

    static bool combined_region(
        mio::smm::SimulationSet<num_regions, mio::osir::InfectionState, closure_order>& stochastic_model,
        mio::smm_moments::Simulation<num_regions, closure_order>& deterministic_model, bool stochastic_used,
        size_t region)
    {
        if (stochastic_used) {
            auto relations = current_smm_relations(stochastic_model, 0.);
            auto current_means = current_smm_means(stochastic_model, 0.);
            double R0          = m_config->lambdas[region] / m_config->gamma *
                        current_means[region * (int)mio::osir::InfectionState::Count +
                                      (int)mio::osir::InfectionState::Susceptible];
            if (((relations[region * (int)mio::osir::InfectionState::Count + (int)mio::osir::InfectionState::Infected] >
                0) && (relations[region * (int)mio::osir::InfectionState::Count + (int)mio::osir::InfectionState::Infected] <
                m_mean_stddev_relation)) || (R0 < 1.0)) {
                return true;
            }
            return false;
        }
        else {
            auto relations = current_moment_relations(deterministic_model, 0.);
            auto current_means = current_moment_means(deterministic_model, 0.);
            double R0          = m_config->lambdas[region] / m_config->gamma *
                        current_means[region * (int)mio::osir::InfectionState::Count +
                                      (int)mio::osir::InfectionState::Susceptible];
            if (((relations[region * (int)mio::osir::InfectionState::Count + (int)mio::osir::InfectionState::Infected] <=
                0) || (relations[region * (int)mio::osir::InfectionState::Count + (int)mio::osir::InfectionState::Infected] >
                m_mean_stddev_relation)) && (R0 > 1.0)) {
                return true;
            }
            return false;
        }
    }

    static bool
    pure_ode(mio::smm::SimulationSet<num_regions, mio::osir::InfectionState, closure_order>& /*stochastic_model*/,
             mio::smm_moments::Simulation<num_regions, closure_order>& /*deterministic_model*/, bool stochastic_used,
             size_t /*region*/)
    {
        if (stochastic_used) {
            return true;
        }
        else {
            return false;
        }
    }

    static bool pure_stochastic(
        mio::smm::SimulationSet<num_regions, mio::osir::InfectionState, closure_order>& /*stochastic_model*/,
        mio::smm_moments::Simulation<num_regions, closure_order>& /*deterministic_model*/, bool stochastic_used,
        size_t /*region*/)
    {
        if (stochastic_used) {
            return false;
        }
        else {
            return true;
        }
    }

    static std::vector<double>
    current_smm_means(mio::smm::SimulationSet<num_regions, mio::osir::InfectionState, closure_order>& sim, double /*t*/)
    {
        return sim.get_last_means();
    }

    static std::vector<double> current_moment_means(mio::smm_moments::Simulation<num_regions, closure_order>& sim,
                                                    double /*t*/)
    {
        return sim.get_last_means();
    }

    static std::vector<double>
    current_smm_var_gradients(mio::smm::SimulationSet<num_regions, mio::osir::InfectionState, closure_order>& sim,
                              double /*t*/)
    {
        return sim.get_last_var_gradients();
    }

    static std::vector<double>
    current_moment_var_gradients(mio::smm_moments::Simulation<num_regions, closure_order>& sim, double /*t*/)
    {
        return sim.get_last_var_gradients();
    }

    static std::vector<double>
    current_smm_mean_gradients(mio::smm::SimulationSet<num_regions, mio::osir::InfectionState, closure_order>& sim,
                               double /*t*/)
    {
        return sim.get_last_mean_gradients();
    }

    static std::vector<double>
    current_moment_mean_gradients(mio::smm_moments::Simulation<num_regions, closure_order>& sim, double /*t*/)
    {
        return sim.get_last_mean_gradients();
    }

    static std::vector<double>
    current_smm_relations(mio::smm::SimulationSet<num_regions, mio::osir::InfectionState, closure_order>& sim,
                          double /*t*/)
    {
        auto last_means = sim.get_last_means();
        auto last_vars  = sim.get_last_vars();
        std::vector<double> relation(last_means.size());
        for (size_t i = 0; i < relation.size(); ++i) {
            relation[i] = last_means[i] > 0 ? std::sqrt(last_vars[i]) / last_means[i] : 0;
        }
        return relation;
    }

    static std::vector<double> current_moment_relations(mio::smm_moments::Simulation<num_regions, closure_order>& sim,
                                                        double /*t*/)
    {
        auto last_means = sim.get_last_means();
        auto last_vars  = sim.get_last_vars();
        std::vector<double> relation(last_means.size());
        for (size_t i = 0; i < relation.size(); ++i) {
            relation[i] = last_means[i] > 0 ? std::sqrt(last_vars[i]) / last_means[i] : 0;
        }
        return relation;
    }
};
} // namespace hybrid

} // namespace mio
#endif // CONDITION_FUNCTIONS_H
