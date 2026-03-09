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
#include "ode_sir/infection_state.h"
#include "simulations/hybrid_simulations/sir_metapop/config/config.cpp"
#include "memilio/utils/logging.h"
#include "simulations/hybrid_simulations/sir_metapop/library/moment_helper.h"
#include "simulations/hybrid_simulations/sir_metapop/library/smm_helper.h"
#include "smm/simulation_set.h"
#include "smm_moments/simulation.h"
#include "models/hybrid/conversion_functions.cpp"
#include "smm_moments/closure_functions.h"
#include "smm/model.h"
#include <cstddef>
#include <limits>
#include <utility>

template <class ModelFrom, class ModelTo>
void exchange_agents(ModelFrom&, ModelTo&, size_t, size_t) = delete;

template <>
void exchange_agents(mio::smm::SimulationSet<2, mio::osir::InfectionState, 3>& model_from,
                     mio::smm::SimulationSet<2, mio::osir::InfectionState, 3>& model_to, size_t /*region_from*/,
                     size_t region_to)
{
    auto& mean_from = model_from.get_mean();
    bool exchange   = false;

    for (size_t comp = 0; comp < static_cast<size_t>(mio::osir::InfectionState::Count); ++comp) {
        double value =
            mean_from.get_last_value()[static_cast<size_t>(mio::osir::InfectionState::Count) * region_to + comp];
        if (value >= 1) {
            if (!exchange) {
                exchange = true;
                break;
            }
        }
    }
    if (!exchange) {
        return;
    }

    // Update simulation results in both regions
    auto& simulations_from = model_from.get_simulations();
    auto& simulations_to   = model_to.get_simulations();
    auto& results_from     = model_from.get_result();
    auto& results_to       = model_to.get_result();
    for (size_t sim = 0; sim < simulations_from.size(); ++sim) {
        auto& to_sim_result    = simulations_to[sim].get_result();
        auto& from_sim_result  = simulations_from[sim].get_result();
        auto& to_sim_result1   = simulations_to[sim].get_result1();
        auto& from_sim_result1 = simulations_from[sim].get_result1();
        auto& to_pop           = simulations_to[sim].get_model().populations;
        auto& from_pop         = simulations_from[sim].get_model().populations;
        for (size_t comp = 0; comp < static_cast<size_t>(mio::osir::InfectionState::Count); ++comp) {
            to_sim_result.get_last_value()[static_cast<size_t>(mio::osir::InfectionState::Count) * region_to + comp] +=
                from_sim_result
                    .get_last_value()[static_cast<size_t>(mio::osir::InfectionState::Count) * region_to + comp];
            from_sim_result.get_last_value()[static_cast<size_t>(mio::osir::InfectionState::Count) * region_to + comp] =
                0;
            to_sim_result1.get_last_value()[static_cast<size_t>(mio::osir::InfectionState::Count) * region_to + comp] +=
                from_sim_result
                    .get_last_value()[static_cast<size_t>(mio::osir::InfectionState::Count) * region_to + comp];
            from_sim_result1
                .get_last_value()[static_cast<size_t>(mio::osir::InfectionState::Count) * region_to + comp] = 0;
            to_pop[{mio::regions::Region(region_to), mio::osecir::InfectionState(comp)}] +=
                from_pop[{mio::regions::Region(region_to), mio::osecir::InfectionState(comp)}];
            from_pop[{mio::regions::Region(region_to), mio::osecir::InfectionState(comp)}] = 0;
            results_to[sim]
                .get_last_value()[static_cast<size_t>(mio::osir::InfectionState::Count) * region_to + comp] +=
                results_from[sim]
                    .get_last_value()[static_cast<size_t>(mio::osir::InfectionState::Count) * region_to + comp];
            results_from[sim]
                .get_last_value()[static_cast<size_t>(mio::osir::InfectionState::Count) * region_to + comp] = 0;
        }
    }

    //Update mean timeseries in both regions
    mean_from.remove_last_time_point();
    model_from.recalculate_last_mean();
    auto& mean_to = model_to.get_mean();
    mean_to.remove_last_time_point();
    model_to.recalculate_last_mean();

    //Update moments in both regions
    auto& moments_from = model_from.get_moments();
    moments_from.remove_last_time_point();
    model_from.recalculate_last_moments();
    auto& moments_to = model_to.get_moments();
    moments_to.remove_last_time_point();
    model_to.recalculate_last_moments();
}

template <>
void exchange_agents(mio::smm_moments::Simulation<2, 3>& model_from, mio::smm_moments::Simulation<2, 3>& model_to,
                     size_t /*region_from*/, size_t region_to)
{
    const int closure_order = 3;
    auto& model_from_result = model_from.get_result();
    auto& model_from_model  = model_from.get_model();
    bool exchange           = false;
    for (size_t comp = 0; comp < static_cast<size_t>(mio::osir::InfectionState::Count); ++comp) {
        size_t comp_index_region_to = model_from_model.populations.get_flat_index(
            {mio::regions::Region(region_to), mio::osir::InfectionState(comp)});
        double value = model_from_result.get_last_value()[comp_index_region_to];
        if (value > 0) {
            if (!exchange) {
                exchange = true;
                break;
            }
        }
    }
    if (!exchange) {
        return;
    }

    auto& model_to_result = model_to.get_result();
    auto& model_to_model  = model_to.get_model();

    // Exchange means
    for (size_t comp = 0; comp < static_cast<size_t>(mio::osir::InfectionState::Count); ++comp) {
        size_t comp_index_region_to = model_from_model.populations.get_flat_index(
            {mio::regions::Region(region_to), mio::osir::InfectionState(comp)});
        model_to_result.get_last_value()[comp_index_region_to] +=
            model_from_result.get_last_value()[comp_index_region_to];
        model_from_result.get_last_value()[comp_index_region_to] = 0;
        model_to_model.populations[{mio::regions::Region(region_to), mio::osecir::InfectionState(comp)}] +=
            model_from_result.get_last_value()[comp_index_region_to];
        model_from_model.populations[{mio::regions::Region(region_to), mio::osecir::InfectionState(comp)}] = 0;
    }

    //Exchange moments
    for (size_t index = 0;
         index < model_to_result.get_last_value().size() - model_to_model.populations.get_num_compartments(); ++index) {
        auto multiindex = model_to_model.moments.unflatten_index(index);
        int order       = std::accumulate(multiindex.begin(), multiindex.end(), 0.0);
        if (order < closure_order) {
            int order_region_to =
                std::accumulate(multiindex.begin() + static_cast<size_t>(mio::osir::InfectionState::Count) * region_to,
                                multiindex.begin() + static_cast<size_t>(mio::osir::InfectionState::Count) * region_to +
                                    static_cast<size_t>(mio::osir::InfectionState::Count),
                                0.0);
            if (order_region_to == 0) {
                continue;
            }
            else if (order == order_region_to) {
                model_to_result.get_last_value()[index] += model_from_result.get_last_value()[index];
                model_from_result.get_last_value()[index] = 0;
            }
            else {
                model_to_result.get_last_value()[index]   = 0;
                model_from_result.get_last_value()[index] = 0;
            }
        }
    }
}

template <>
void exchange_agents(mio::smm::SimulationSet<2, mio::osir::InfectionState, 3>& model_from,
                     mio::smm_moments::Simulation<2, 3>& model_to, size_t /*region_from*/, size_t region_to)
{
    auto& mean_smm_set = model_from.get_mean();
    bool exchange      = false;

    for (size_t comp = 0; comp < static_cast<size_t>(mio::osir::InfectionState::Count); ++comp) {
        double value =
            mean_smm_set.get_last_value()[static_cast<size_t>(mio::osir::InfectionState::Count) * region_to + comp];
        if (value > 0) {
            if (!exchange) {
                exchange = true;
                break;
            }
        }
    }
    if (!exchange) {
        return;
    }

    auto& moment_result = model_to.get_result();
    auto& moment_model  = model_to.get_model();

    // Update means of moment model
    for (size_t comp = 0; comp < static_cast<size_t>(mio::osir::InfectionState::Count); ++comp) {
        size_t comp_index_region_to =
            moment_model.populations.get_flat_index({mio::regions::Region(region_to), mio::osir::InfectionState(comp)});
        moment_result.get_last_value()[comp_index_region_to] += mean_smm_set.get_last_value()[comp_index_region_to];
        moment_model.populations[{mio::regions::Region(region_to), mio::osecir::InfectionState(comp)}] +=
            mean_smm_set.get_last_value()[comp_index_region_to];
    }

    //Update moments of moment model
    auto moment_indices_smm_set = model_from.get_moment_indices();
    auto& moments_smm_set       = model_from.get_moments();
    for (size_t i = 0; i < moment_indices_smm_set.size(); ++i) {
        int order           = std::accumulate(moment_indices_smm_set[i].begin(), moment_indices_smm_set[i].end(), 0.0);
        int order_region_to = std::accumulate(
            moment_indices_smm_set[i].begin() + static_cast<size_t>(mio::osir::InfectionState::Count) * region_to,
            moment_indices_smm_set[i].begin() + static_cast<size_t>(mio::osir::InfectionState::Count) * region_to +
                static_cast<size_t>(mio::osir::InfectionState::Count),
            0.0);
        if (order_region_to != 0) {
            if (order_region_to == order) {
                auto index = moment_model.moments.flatten_index(moment_indices_smm_set[i]);
                moment_result.get_last_value()[moment_model.populations.get_num_compartments() + index] +=
                    moments_smm_set.get_last_value()[i];
            }
        }
    }

    auto& smm_set_sims    = model_from.get_simulations();
    auto& smm_set_results = model_from.get_result();
    //Reset results and simulations of simulation set
    for (size_t sim = 0; sim < smm_set_sims.size(); ++sim) {
        auto& sim_result = smm_set_sims[sim].get_result();
        auto& sim_pop    = smm_set_sims[sim].get_model().populations;
        for (size_t comp = 0; comp < static_cast<size_t>(mio::osir::InfectionState::Count); ++comp) {
            sim_result.get_last_value()[static_cast<size_t>(mio::osir::InfectionState::Count) * region_to + comp] = 0;
            sim_pop[{mio::regions::Region(region_to), mio::osecir::InfectionState(comp)}]                         = 0;
            smm_set_results[sim]
                .get_last_value()[static_cast<size_t>(mio::osir::InfectionState::Count) * region_to + comp] = 0;
        }
    }

    //Reset and update mean and moment ts of simulation set
    mean_smm_set.remove_last_time_point();
    model_from.recalculate_last_mean();
    moments_smm_set.remove_last_time_point();
    model_from.recalculate_last_moments();
}

template <>
void exchange_agents(mio::smm_moments::Simulation<2, 3>& model_from,
                     mio::smm::SimulationSet<2, mio::osir::InfectionState, 3>& model_to, size_t /*region_from*/,
                     size_t region_to)
{
    const size_t num_regions = 2;
    auto& moment_result      = model_from.get_result();
    auto& moment_model       = model_from.get_model();
    bool exchange            = false;
    for (size_t comp = 0; comp < static_cast<size_t>(mio::osir::InfectionState::Count); ++comp) {
        size_t comp_index_region_to =
            moment_model.populations.get_flat_index({mio::regions::Region(region_to), mio::osir::InfectionState(comp)});
        double value = moment_result.get_last_value()[comp_index_region_to];
        if (value >= 1) {
            if (!exchange) {
                exchange = true;
                break;
            }
        }
    }
    if (!exchange) {
        return;
    }
    // For each smm simulation draw number of individuals coming to "to" region from mean and variances of moments model for every compartment
    auto& smm_set_sims    = model_to.get_simulations();
    auto& smm_set_results = model_to.get_result();
    for (size_t sim = 0; sim < smm_set_sims.size(); ++sim) {
        auto& sim_result  = smm_set_sims[sim].get_result();
        auto& sim_result1 = smm_set_sims[sim].get_result1();
        auto& sim_pop     = smm_set_sims[sim].get_model().populations;
        for (size_t comp = 0; comp < static_cast<size_t>(mio::osir::InfectionState::Count); ++comp) {
            double mean =
                moment_result
                    .get_last_value()[region_to * static_cast<size_t>(mio::osir::InfectionState::Count) + comp];
            std::array<int, static_cast<size_t>(mio::osir::InfectionState::Count) * num_regions> indices_var;
            indices_var.fill(0);
            indices_var[region_to * static_cast<size_t>(mio::osir::InfectionState::Count) + comp] = 2;
            double var       = moment_result.get_last_value()[moment_model.populations.get_num_compartments() +
                                                        moment_model.moments.flatten_index(indices_var)];
            double sim_value = std::max(0., std::round(mio::NormalDistribution<double>::get_instance()(
                                                smm_set_sims[sim].get_model().get_rng(), mean, var)));
            sim_result.get_last_value()[static_cast<size_t>(mio::osir::InfectionState::Count) * region_to + comp] +=
                sim_value;
            sim_result1.get_last_value()[static_cast<size_t>(mio::osir::InfectionState::Count) * region_to + comp] +=
                sim_value;
            sim_pop[{mio::regions::Region(region_to), mio::osecir::InfectionState(comp)}] += sim_value;
            smm_set_results[sim]
                .get_last_value()[static_cast<size_t>(mio::osir::InfectionState::Count) * region_to + comp] +=
                sim_value;
        }
    }

    //Reset and update mean and moment ts of simulation set
    model_to.get_mean().remove_last_time_point();
    model_to.recalculate_last_mean();
    model_to.get_moments().remove_last_time_point();
    model_to.recalculate_last_moments();

    // Set means and moments in region to in moment model to zero
    for (size_t index = 0; index < static_cast<size_t>(moment_result.get_num_elements()); ++index) {
        if (index >= region_to * static_cast<size_t>(mio::osir::InfectionState::Count) &&
            index < region_to * static_cast<size_t>(mio::osir::InfectionState::Count) +
                        static_cast<size_t>(mio::osir::InfectionState::Count)) {
            moment_result.get_last_value()[index] = 0;
        }
        if (index >= moment_model.populations.get_num_compartments()) {
            auto multiindex =
                moment_model.moments.unflatten_index(index - moment_model.populations.get_num_compartments());
            int order_region_to =
                std::accumulate(multiindex.begin() + static_cast<size_t>(mio::osir::InfectionState::Count) * region_to,
                                multiindex.begin() + static_cast<size_t>(mio::osir::InfectionState::Count) * region_to +
                                    static_cast<size_t>(mio::osir::InfectionState::Count),
                                0.0);
            if (order_region_to > 0) {
                moment_result.get_last_value()[index] = 0;
            }
        }
    }
}

template <size_t num_regions, size_t closure_order>
class SpatialHybridSim
{

public:
    using SMMSetSim  = mio::smm::SimulationSet<num_regions, mio::osir::InfectionState, closure_order>;
    using MomentsSim = mio::smm_moments::Simulation<num_regions, closure_order>;
    using TemporalHybridSim =
        mio::hybrid::TemporalHybridSimulation<SMMSetSim, MomentsSim,
                                              std::pair<std::vector<double>, std::vector<double>>,
                                              std::pair<std::vector<double>, std::vector<double>>>;
    SpatialHybridSim(const Config::Config& config, size_t num_runs, double min_step_size, double dt_switch, double t0)
        : m_simulations()
        , m_t(t0)
        , m_dt(dt_switch)
    {
        // Initialize one model per region and delete all adoption rates that are not in the region and all transition rates that are not outgoing the region
        m_simulations.reserve(num_regions);
        for (size_t region = 0; region < num_regions; ++region) {
            m_simulations.emplace_back(
                initialize_temporal_hybrid_model(config, num_runs, min_step_size, dt_switch, region));
        }
    }

    void advance(double tmax)
    {
        //Condition: we don't switch if the relation of mean and stddev is bigger than 60%
        const auto condition_func = [](std::pair<std::vector<double>, std::vector<double>>& result_smm,
                                       std::pair<std::vector<double>, std::vector<double>>& /*result_moments*/,
                                       bool smm_used) {
            if (smm_used) {
                for (size_t r = 0; r < num_regions; ++r) {
                    auto var_infected_gradient = result_smm.first[r * (int)mio::osir::InfectionState::Count +
                                                                  (int)mio::osir::InfectionState::Infected];
                    auto relation_infected     = result_smm.second[r * (int)mio::osir::InfectionState::Count +
                                                               (int)mio::osir::InfectionState::Infected];
                    if (var_infected_gradient < -1 && relation_infected < 0.6 && relation_infected > 0) {
                        return true;
                    }
                }
            }

            return false;
        };
        while (m_t <= tmax) {
            for (auto& region_sim : m_simulations) {
                region_sim.advance(m_t + m_dt, condition_func);
            }
            for (size_t region_from = 0; region_from < num_regions; ++region_from) {
                for (size_t region_to = 0; region_to < num_regions; ++region_to) {
                    if (region_from != region_to) {
                        exchange(m_simulations[region_from], region_from, m_simulations[region_to], region_to);
                        auto a = m_simulations[region_from].get_model2().get_result();
                        mio::unused(a);
                    }
                }
            }
            m_t += m_dt;
        }
    }

    const TemporalHybridSim& get_sim_by_region(size_t region) const
    {
        return m_simulations[region];
    }
    TemporalHybridSim& get_sim_by_region(size_t region)
    {
        return m_simulations[region];
    }

private:
    void exchange(TemporalHybridSim& model_from, size_t region_from, TemporalHybridSim& model_to, size_t region_to)
    {
        if (model_from.using_model1()) {
            auto& model_from_used = model_from.get_model1();
            if (model_to.using_model1()) {
                auto& model_to_used = model_to.get_model1();
                exchange_agents(model_from_used, model_to_used, region_from, region_to);
            }
            else {
                auto& model_to_used = model_to.get_model2();
                exchange_agents(model_from_used, model_to_used, region_from, region_to);
            }
        }
        else {
            auto& model_from_used = model_from.get_model2();
            if (model_to.using_model1()) {
                auto& model_to_used = model_to.get_model1();
                exchange_agents(model_from_used, model_to_used, region_from, region_to);
            }
            else {
                auto& model_to_used = model_to.get_model2();
                exchange_agents(model_from_used, model_to_used, region_from, region_to);
            }
        }
    }

    mio::smm::Model<ScalarType, num_regions, mio::osir::InfectionState>
    initialize_smm_for_region(const Config::Config& config, size_t region)
    {
        assert(config.num_regions == num_regions);
        mio::smm::Model<ScalarType, num_regions, mio::osir::InfectionState> model;

        // Initialize populations
        for (size_t r = 0; r < config.num_regions; ++r) {
            model.populations[{mio::regions::Region(r), mio::osir::InfectionState::Susceptible}] = 0;
            model.populations[{mio::regions::Region(r), mio::osir::InfectionState::Infected}]    = 0;
            model.populations[{mio::regions::Region(r), mio::osir::InfectionState::Recovered}]   = 0;

            if (r == region) {
                model.populations[{mio::regions::Region(r), mio::osir::InfectionState::Susceptible}] =
                    config.total_populations[r];
            }
        }
        // Set initially infected
        for (size_t i = 0; i < config.I0s.size(); ++i) {
            int region_id = config.I0s[i].first;
            double I0     = config.I0s[i].second;
            if (static_cast<size_t>(region_id) == region) {
                model.populations[{mio::regions::Region(region_id), mio::osir::InfectionState::Infected}] = I0;
                model.populations[{mio::regions::Region(region_id), mio::osir::InfectionState::Susceptible}] =
                    config.total_populations[region_id] -
                    model.populations[{mio::regions::Region(region_id), mio::osir::InfectionState::Infected}] -
                    model.populations[{mio::regions::Region(region_id), mio::osir::InfectionState::Recovered}];
                break;
            }
        }

        std::vector<mio::AdoptionRate<ScalarType, mio::osir::InfectionState>> adoption_rates;
        for (size_t r = 0; r < config.num_regions; ++r) {
            // Second-order adoption rate lambda is region dependent
            adoption_rates.push_back({mio::osir::InfectionState::Susceptible,
                                      mio::osir::InfectionState::Infected,
                                      mio::regions::Region(r),
                                      config.lambdas[r],
                                      {{mio::osir::InfectionState::Infected, 1.}}});
            // Recovery rate gamma is the same for all regions
            adoption_rates.push_back({mio::osir::InfectionState::Infected,
                                      mio::osir::InfectionState::Recovered,
                                      mio::regions::Region(r),
                                      config.gamma,
                                      {}});
        }

        model.parameters.template get<mio::smm::AdoptionRates<ScalarType, mio::osir::InfectionState>>() =
            adoption_rates;

        std::vector<mio::smm::TransitionRate<ScalarType, mio::osir::InfectionState>> transition_rates;
        for (const auto& rate : config.transition_rates) {
            if (rate.from == mio::regions::Region(region)) {
                transition_rates.push_back(rate);
            }
        }
        model.parameters.template get<mio::smm::TransitionRates<ScalarType, mio::osir::InfectionState>>() =
            transition_rates;

        return model;
    }

    mio::smm_moments::Model<num_regions, closure_order> initialize_moments_for_region(
        Eigen::Array<double, Eigen::Dynamic, 1>& expected_values_init,
        Eigen::Array<double, Eigen::Dynamic, 1>& moments_init, const Config::Config& config,
        typename mio::smm_moments::Model<num_regions, closure_order>::ClosureFunctionType closure_func, size_t region)
    {
        mio::smm_moments::Model<num_regions, closure_order> model(closure_func);
        // Check whether initial expected values and moments have the correct size
        assert(expected_values_init.rows() == num_regions * static_cast<size_t>(mio::osir::InfectionState::Count) &&
               "Initial expected values do not have correct size");
        assert(moments_init.rows() == model.moments.moments().rows() && "Initial moments do not have correct size");

        // Set spatial transition rates
        for (auto& rate : config.transition_rates) {
            if (rate.from == mio::regions::Region(region)) {
                model.parameters.template get<mio::smm_moments::TransitionRate>()[{rate.status, rate.from, rate.to}] =
                    rate.factor;
            }
        }

        for (size_t r = 0; r < config.num_regions; ++r) {
            // Set recovery rate
            model.parameters.template get<mio::smm_moments::RecoveryRate>()[mio::regions::Region(r)] = config.gamma;
            // Set transmission rates
            model.parameters.template get<mio::smm_moments::TransmissionRate>()[mio::regions::Region(r)] =
                config.lambdas[r];
        }

        for (size_t r = 0; r < num_regions; ++r) {
            // Set initial expected values
            model.populations[{mio::regions::Region(r), mio::osir::InfectionState::Susceptible}] = 0;
            model.populations[{mio::regions::Region(r), mio::osir::InfectionState::Infected}]    = 0;
            model.populations[{mio::regions::Region(r), mio::osir::InfectionState::Recovered}]   = 0;
            if (r == region) {
                model.populations[{mio::regions::Region(r), mio::osir::InfectionState::Susceptible}] =
                    expected_values_init[r * static_cast<size_t>(mio::osir::InfectionState::Count) +
                                         static_cast<size_t>(mio::osir::InfectionState::Susceptible)];
                model.populations[{mio::regions::Region(r), mio::osir::InfectionState::Infected}] =
                    expected_values_init[r * static_cast<size_t>(mio::osir::InfectionState::Count) +
                                         static_cast<size_t>(mio::osir::InfectionState::Infected)];
                model.populations[{mio::regions::Region(r), mio::osir::InfectionState::Recovered}] =
                    expected_values_init[r * static_cast<size_t>(mio::osir::InfectionState::Count) +
                                         static_cast<size_t>(mio::osir::InfectionState::Recovered)];
            }
        }
        // Set initial moments
        model.moments.moments() = moments_init;
        return model;
    }

    TemporalHybridSim initialize_temporal_hybrid_model(const Config::Config& config, size_t num_runs,
                                                       double min_step_size, double dt_switch, size_t region)
    {
        // Initialize smm
        auto smm_model = initialize_smm_for_region(config, region);
        // Initialize smm simulation set
        auto sim_set = SMMSetSim(num_runs, smm_model, config.t0, config.dt);
        // Initialize moment model
        // Initial expected values
        Eigen::Array<double, Eigen::Dynamic, 1> expected_values_init(
            num_regions * static_cast<size_t>(mio::osir::InfectionState::Count));
        expected_values_init.setZero();
        for (size_t r = 0; r < num_regions; ++r) {
            for (size_t s = 0; s < static_cast<size_t>(mio::osir::InfectionState::Count); ++s) {
                expected_values_init[r * static_cast<size_t>(mio::osir::InfectionState::Count) + s] =
                    smm_model.populations[{mio::regions::Region(r), mio::osir::InfectionState(s)}];
            }
        }
        auto closure_func = &mio::smm_moments::truncation_closure<num_regions, closure_order>;
        // Initial moments are all zero
        MomentArray<static_cast<size_t>(mio::osir::InfectionState::Count), num_regions, closure_order> moments_array;
        std::array<int, static_cast<size_t>(mio::osir::InfectionState::Count) * num_regions> zero_index;
        zero_index.fill(0);
        moments_array.moments()[moments_array.flatten_index(zero_index)] = 1.0;
        auto moment_model =
            initialize_moments_for_region(expected_values_init, moments_array.moments(), config, closure_func, region);

        // Initialize moment simulation
        auto sim_moments = MomentsSim(moment_model, config.t0, config.dt);
        // Set maximum dt of integrator to interpolation time points
        sim_moments.get_integrator_core().get_dt_max() = config.dt;
        if (min_step_size > 0) {
            sim_moments.get_integrator_core().get_dt_min() = min_step_size;
        }

        // Define result function smm
        const auto result_fct_smm =
            [](mio::smm::SimulationSet<num_regions, mio::osir::InfectionState, closure_order>& sim, double /*t*/) {
                auto last_var_gradients = sim.get_last_var_gradients();
                auto last_means         = sim.get_last_means();
                auto last_vars          = sim.get_last_vars();
                std::vector<double> relation(last_means.size());
                for (size_t i = 0; i < relation.size(); ++i) {
                    relation[i] = last_means[i] > 0 ? std::sqrt(last_vars[i]) / last_means[i] : 0;
                }
                return std::pair<std::vector<double>, std::vector<double>>{last_var_gradients, relation};
            };

        const auto result_fct_moments = [](mio::smm_moments::Simulation<num_regions, closure_order>& sim,
                                           double /*t*/) {
            auto last_var_gradients = sim.get_last_var_gradients();
            auto last_means         = sim.get_last_means();
            auto last_vars          = sim.get_last_vars();
            std::vector<double> relation(last_means.size());
            for (size_t i = 0; i < relation.size(); ++i) {
                relation[i] = last_means[i] > 0 ? std::sqrt(last_vars[i]) / last_means[i] : 0;
            }
            return std::pair<std::vector<double>, std::vector<double>>{last_var_gradients, relation};
        };

        return TemporalHybridSim(std::move(sim_set), std::move(sim_moments), result_fct_smm, result_fct_moments, true,
                                 config.t0, dt_switch);
    }

    std::vector<TemporalHybridSim> m_simulations;
    double m_t;
    double m_dt;
};

int main()
{
    mio::set_log_level(mio::LogLevel::warn);
    const size_t num_runs      = 10000;
    double dt_switch           = 1.;
    const size_t closure_order = 3;
    const auto config          = Config::get_config(Config::ConfigType::Config2regionsk1);
    const size_t num_regions   = 2;
    double min_step_size       = 0.0001;

    if (num_regions != config.num_regions) {
        mio::log_error("Number of regions doesn't match number of regions in config.");
    }

    std::string save_file = Config::SAVE_DIR + "Spatial-Hybrid2/";
    save_file += config.name;

    auto created_directory = mio::create_directory(save_file);
    if (!created_directory) {
        printf("%s\n", created_directory.error().formatted_message().c_str());
        return -1;
    }

    save_file += "/" + Config::switch_condition_string[1];
    created_directory = mio::create_directory(save_file);
    if (!created_directory) {
        printf("%s\n", created_directory.error().formatted_message().c_str());
        return -1;
    }

    save_file += "/" + Config::closure_string[0];
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

    // Second temp hybrid model for second region and delete the rates per model in every region
    SpatialHybridSim<num_regions, closure_order> spatial_hybrid_sim(config, num_runs, min_step_size, dt_switch,
                                                                    config.t0);

    mio::timing::BasicTimer timer;
    timer.start();
    spatial_hybrid_sim.advance(config.tmax);
    timer.stop();

    //Calculate outputs per region
    for (size_t r = 0; r < num_regions; ++r) {
        auto region_sim = spatial_hybrid_sim.get_sim_by_region(r);
        // Calculate moments and expected values of sim set
        auto means_smm       = region_sim.get_model1().get_mean();
        auto moments_smm     = region_sim.get_model1().get_moments();
        auto done            = means_smm.export_csv(save_file + std::to_string(r) + "_means_smm.csv");
        done                 = moments_smm.export_csv(save_file + std::to_string(r) + "_moments_smm.csv",
                                                      region_sim.get_model1().get_moment_names());
        auto expected_values = region_sim.get_model2().get_expected_values_time_series();
        auto moments_moments = region_sim.get_model2().get_moment_time_series(closure_order);
        done = expected_values.export_csv(save_file + std::to_string(r) + "_expected_values_moments.csv");
        done = moments_moments.first.export_csv(save_file + std::to_string(r) + "_moments_moments.csv",
                                                moments_moments.second);

        // Get percentiles of smm simulation set
        std::vector<std::vector<mio::TimeSeries<double>>> sim_results;
        auto& all_results = region_sim.get_model1().get_result();
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
        auto finished = p05[0].export_csv(save_file + std::to_string(r) + "_smm_p05.csv");
        finished      = p25[0].export_csv(save_file + std::to_string(r) + "_smm_p25.csv");
        finished      = p50[0].export_csv(save_file + std::to_string(r) + "_smm_p50.csv");
        finished      = p75[0].export_csv(save_file + std::to_string(r) + "_smm_p75.csv");

        // Save merged time series
        auto hybrid_result_means   = mio::merge_time_series(means_smm, expected_values).value();
        auto hybrid_result_moments = mio::merge_time_series(moments_smm, moments_moments.first).value();

        int num_steps = static_cast<int>(config.tmax / config.dt) + 1;
        std::vector<double> interpolation_tps(num_steps);

        for (int i = 0; i < num_steps; ++i) {
            interpolation_tps[i] = i * config.dt;
        }
        hybrid_result_means   = mio::interpolate_simulation_result(hybrid_result_means, interpolation_tps);
        hybrid_result_moments = mio::interpolate_simulation_result(hybrid_result_moments, interpolation_tps);
        done                  = hybrid_result_means.export_csv(save_file + std::to_string(r) + "_means.csv");
        done = hybrid_result_moments.export_csv(save_file + std::to_string(r) + "_moments.csv", moments_moments.second);
    }

    mio::TimeSeries<double> total_time(1);
    Eigen::VectorXd time = Eigen::VectorXd::Constant(1, mio::timing::time_in_seconds(timer.get_elapsed_time()));
    total_time.add_time_point(0., time);
    auto finished_time = total_time.export_csv(save_file + "total_time.csv", {"Runtime"});

    std::cout << "Spatial-hybrid Elapsed time: " << timer.get_elapsed_time() << std::endl << std::flush;

    return 0;
}
