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
#include "smm/simulation_set.h"
#include "smm_moments/simulation.h"
#include "smm_moments/closure_functions.h"
#include "smm/model.h"
#include "simulations/hybrid_simulations/sir_metapop/config/config.h"
#include "hybrid/conversion_functions.h"
#include "hybrid/exchange_functions.h"

namespace mio
{
namespace hybrid
{

template <size_t num_regions, size_t closure_order>
class SpatialHybridSimulation
{

public:
    using SMMSetSim = mio::smm::SimulationSet<num_regions, mio::osir::InfectionState, closure_order>;
    using MomentSim = mio::smm_moments::Simulation<num_regions, closure_order>;
    using TemporalHybridSim =
        mio::hybrid::TemporalHybridSimulation<SMMSetSim, MomentSim, std::pair<std::vector<double>, std::vector<double>>,
                                              std::pair<std::vector<double>, std::vector<double>>>;
    SpatialHybridSimulation(const Config::Config& config, size_t num_runs, double min_step_size, double dt_switch,
                            double t0)
        : m_simulations()
        , m_t(t0)
        , m_dt(dt_switch)
    {
        // Initialize one model per region and delete all transition rates that are not outgoing the region
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
        auto sim_moments = MomentSim(moment_model, config.t0, config.dt);
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

} // namespace hybrid

} // namespace mio

#endif //MIO_SPATIAL_HYBRID_MODEL_H
