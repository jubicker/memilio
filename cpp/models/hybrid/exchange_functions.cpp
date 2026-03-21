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

#include "hybrid/exchange_functions.h"

namespace mio
{
namespace hybrid
{

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
            to_pop[{mio::regions::Region(region_to), mio::osir::InfectionState(comp)}] +=
                from_pop[{mio::regions::Region(region_to), mio::osir::InfectionState(comp)}];
            from_pop[{mio::regions::Region(region_to), mio::osir::InfectionState(comp)}] = 0;
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
        model_to_model.populations[{mio::regions::Region(region_to), mio::osir::InfectionState(comp)}] +=
            model_from_result.get_last_value()[comp_index_region_to];
        model_from_model.populations[{mio::regions::Region(region_to), mio::osir::InfectionState(comp)}] = 0;
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
        moment_model.populations[{mio::regions::Region(region_to), mio::osir::InfectionState(comp)}] +=
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
            sim_pop[{mio::regions::Region(region_to), mio::osir::InfectionState(comp)}]                           = 0;
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
            sim_pop[{mio::regions::Region(region_to), mio::osir::InfectionState(comp)}] += sim_value;
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

} //namespace hybrid

} //namespace mio
