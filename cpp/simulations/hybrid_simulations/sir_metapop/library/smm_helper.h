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
#ifndef SMM_HELPER_H
#define SMM_HELPER_H

#include "hybrid_simulations/sir_metapop/config/config.h"
#include "memilio/utils/time_series.h"
#include "ode_sir/infection_state.h"
#include "smm/model.h"
#include "moment_array.h"
#include <cstddef>
#include <string>
#include <utility>

/**
* @brief Initializes an SMM with given regions and the S-I-R infection states with the given config.
* @param[in] config The configuration containing the model parameters (adoption and transition rates), initial populations and initially infected per region.
*/
template <size_t NumRegions>
mio::smm::Model<NumRegions, mio::osir::InfectionState> initialize_model(const Config::Config& config)
{
    assert(config.num_regions == NumRegions);
    mio::smm::Model<NumRegions, mio::osir::InfectionState> model;

    // Initialize populations
    for (size_t r = 0; r < config.num_regions; ++r) {
        model.populations[{mio::regions::Region(r), mio::osir::InfectionState::Susceptible}] =
            config.total_populations[r];
        model.populations[{mio::regions::Region(r), mio::osir::InfectionState::Infected}]  = 0;
        model.populations[{mio::regions::Region(r), mio::osir::InfectionState::Recovered}] = 0;
    }
    // Set initially infected
    for (size_t i = 0; i < config.I0s.size(); ++i) {
        int region_id = config.I0s[i].first;
        double I0     = config.I0s[i].second;
        model.populations[{mio::regions::Region(region_id), mio::osir::InfectionState::Infected}] = I0;
        model.populations[{mio::regions::Region(region_id), mio::osir::InfectionState::Susceptible}] =
            config.total_populations[region_id] -
            model.populations[{mio::regions::Region(region_id), mio::osir::InfectionState::Infected}] -
            model.populations[{mio::regions::Region(region_id), mio::osir::InfectionState::Recovered}];
    }

    std::vector<mio::AdoptionRate<mio::osir::InfectionState>> adoption_rates;
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
    model.parameters.template get<mio::smm::AdoptionRates<mio::osir::InfectionState>>()   = adoption_rates;
    model.parameters.template get<mio::smm::TransitionRates<mio::osir::InfectionState>>() = config.transition_rates;

    return model;
}

/**
* @brief Calculates the moment defined by indices for the given realizations of a stochastic variable.
* @param[in] values The realizations of the stochastic variable (i.e. S, I, R for each region).
* @param[in] indices The indices defining the moment to be calculated.
*/
template <size_t NumRegions>
double
calculate_moment(const Eigen::Matrix<ScalarType, Eigen::Dynamic,
                                     static_cast<size_t>(mio::osir::InfectionState::Count) * NumRegions>& values,
                 const std::array<int, static_cast<size_t>(mio::osir::InfectionState::Count) * NumRegions>& indices)
{
    // values has 3 columns per region: S, I, R
    Eigen::Matrix<ScalarType, 1, static_cast<size_t>(mio::osir::InfectionState::Count)* NumRegions> means =
        values.colwise().mean();
    double moment = 0.0;
    for (int i = 0; i < values.rows(); ++i) {
        double summand = 1.0;
        for (size_t r = 0; r < NumRegions; ++r) {
            summand *= std::pow(values(i, r * static_cast<size_t>(mio::osir::InfectionState::Count) +
                                              static_cast<size_t>(mio::osir::InfectionState::Susceptible)) -
                                    means(r * static_cast<size_t>(mio::osir::InfectionState::Count) +
                                          static_cast<size_t>(mio::osir::InfectionState::Susceptible)),
                                indices[r * static_cast<size_t>(mio::osir::InfectionState::Count) +
                                        static_cast<size_t>(mio::osir::InfectionState::Susceptible)]) *
                       std::pow(values(i, r * static_cast<size_t>(mio::osir::InfectionState::Count) +
                                              static_cast<size_t>(mio::osir::InfectionState::Infected)) -
                                    means(r * static_cast<size_t>(mio::osir::InfectionState::Count) +
                                          static_cast<size_t>(mio::osir::InfectionState::Infected)),
                                indices[r * static_cast<size_t>(mio::osir::InfectionState::Count) +
                                        static_cast<size_t>(mio::osir::InfectionState::Infected)]) *
                       std::pow(values(i, r * static_cast<size_t>(mio::osir::InfectionState::Count) +
                                              static_cast<size_t>(mio::osir::InfectionState::Recovered)) -
                                    means(r * static_cast<size_t>(mio::osir::InfectionState::Count) +
                                          static_cast<size_t>(mio::osir::InfectionState::Recovered)),
                                indices[r * static_cast<size_t>(mio::osir::InfectionState::Count) +
                                        static_cast<size_t>(mio::osir::InfectionState::Recovered)]);
        }
        moment += summand;
    }
    moment /= static_cast<double>(values.rows());
    return moment;
}

/**
* @brief Calculates the mean timeseries from a vector of time series i.e. at every time point the mean values of all runs at that time point.
* @param[in] sim_result Vector with simulation output (TimeSeries) of an SIR model with NumRegions regions for each run.
* @return Mean timeseries and the corresponding column names.
*/
template <int NumRegions>
std::pair<mio::TimeSeries<ScalarType>, std::vector<std::string>>
calculate_means_from_sim(const std::vector<mio::TimeSeries<ScalarType>>& sim_results)
{
    constexpr int num_states  = static_cast<int>(mio::osir::InfectionState::Count);
    constexpr int num_regions = static_cast<int>(NumRegions);

    Eigen::Matrix<ScalarType, num_states * num_regions, 1> means;
    means.setZero();
    std::vector<std::string> names(static_cast<size_t>(mio::osir::InfectionState::Count) * NumRegions);
    mio::TimeSeries<double> mean_ts(static_cast<size_t>(mio::osir::InfectionState::Count) * NumRegions);
    for (int t = 0; t < sim_results[0].get_num_time_points(); ++t) {
        Eigen::Matrix<ScalarType, Eigen::Dynamic, static_cast<size_t>(mio::osir::InfectionState::Count)* NumRegions>
            result =
                Eigen::Matrix<ScalarType, Eigen::Dynamic,
                              static_cast<size_t>(mio::osir::InfectionState::Count) *
                                  NumRegions>::Zero(sim_results.size(),
                                                    static_cast<size_t>(mio::osir::InfectionState::Count) * NumRegions);
        for (size_t i = 0; i < sim_results.size(); i++) {
            for (size_t r = 0; r < NumRegions; ++r) {
                result(i, static_cast<size_t>(mio::osir::InfectionState::Count) * r + 0) =
                    sim_results[i].get_value(t)[r * static_cast<size_t>(mio::osir::InfectionState::Count) +
                                                static_cast<size_t>(mio::osir::InfectionState::Susceptible)];
                result(i, static_cast<size_t>(mio::osir::InfectionState::Count) * r + 1) =
                    sim_results[i].get_value(t)[r * static_cast<size_t>(mio::osir::InfectionState::Count) +
                                                static_cast<size_t>(mio::osir::InfectionState::Infected)];
                result(i, static_cast<size_t>(mio::osir::InfectionState::Count) * r + 2) =
                    sim_results[i].get_value(t)[r * static_cast<size_t>(mio::osir::InfectionState::Count) +
                                                static_cast<size_t>(mio::osir::InfectionState::Recovered)];
            }
        }
        //Expected values
        for (size_t r = 0; r < NumRegions; ++r) {
            means[r * static_cast<size_t>(mio::osir::InfectionState::Count) +
                  static_cast<size_t>(mio::osir::InfectionState::Susceptible)] =
                result
                    .col(static_cast<size_t>(mio::osir::InfectionState::Count) * r +
                         static_cast<size_t>(mio::osir::InfectionState::Susceptible))
                    .mean(); // muS
            means[r * static_cast<size_t>(mio::osir::InfectionState::Count) +
                  static_cast<size_t>(mio::osir::InfectionState::Infected)] =
                result
                    .col(static_cast<size_t>(mio::osir::InfectionState::Count) * r +
                         static_cast<size_t>(mio::osir::InfectionState::Infected))
                    .mean(); // muI
            means[r * static_cast<size_t>(mio::osir::InfectionState::Count) +
                  static_cast<size_t>(mio::osir::InfectionState::Recovered)] =
                result
                    .col(static_cast<size_t>(mio::osir::InfectionState::Count) * r +
                         static_cast<size_t>(mio::osir::InfectionState::Recovered))
                    .mean(); // muR
        }

        mean_ts.add_time_point(sim_results[0].get_time(t), means);
    }
    for (size_t r = 0; r < NumRegions; ++r) {
        names[r * static_cast<size_t>(mio::osir::InfectionState::Count) +
              static_cast<size_t>(mio::osir::InfectionState::Susceptible)] = "muS_r" + std::to_string(int(r));
        names[r * static_cast<size_t>(mio::osir::InfectionState::Count) +
              static_cast<size_t>(mio::osir::InfectionState::Infected)]    = "muI_r" + std::to_string(int(r));
        names[r * static_cast<size_t>(mio::osir::InfectionState::Count) +
              static_cast<size_t>(mio::osir::InfectionState::Recovered)]   = "muR_r" + std::to_string(int(r));
    }
    return std::make_pair(mean_ts, names);
}

/**
* @brief Calculates a time series for every moment up to order Order from given simulation results.
* @param[in] sim_result Vector with simulation output (TimeSeries) of an SIR model with NumRegions regions for each run.
* @return Timeseries of all moments and the corresponding column names.
*/
template <int Order, int NumRegions>
std::pair<mio::TimeSeries<ScalarType>, std::vector<std::string>>
calculate_moments_from_sim(const std::vector<mio::TimeSeries<ScalarType>>& sim_results)
{
    MomentArray<static_cast<size_t>(mio::osir::InfectionState::Count), NumRegions, Order> moments;
    mio::TimeSeries<double> moment_ts(moments.moments_up_to_order(Order).size());
    for (int t = 0; t < sim_results[0].get_num_time_points(); ++t) {
        Eigen::Matrix<ScalarType, Eigen::Dynamic, static_cast<size_t>(mio::osir::InfectionState::Count)* NumRegions>
            result =
                Eigen::Matrix<ScalarType, Eigen::Dynamic,
                              static_cast<size_t>(mio::osir::InfectionState::Count) *
                                  NumRegions>::Zero(sim_results.size(),
                                                    static_cast<size_t>(mio::osir::InfectionState::Count) * NumRegions);
        for (size_t i = 0; i < sim_results.size(); i++) {
            for (size_t r = 0; r < NumRegions; ++r) {
                result(i, static_cast<size_t>(mio::osir::InfectionState::Count) * r + 0) =
                    sim_results[i].get_value(t)[r * static_cast<size_t>(mio::osir::InfectionState::Count) +
                                                static_cast<size_t>(mio::osir::InfectionState::Susceptible)];
                result(i, static_cast<size_t>(mio::osir::InfectionState::Count) * r + 1) =
                    sim_results[i].get_value(t)[r * static_cast<size_t>(mio::osir::InfectionState::Count) +
                                                static_cast<size_t>(mio::osir::InfectionState::Infected)];
                result(i, static_cast<size_t>(mio::osir::InfectionState::Count) * r + 2) =
                    sim_results[i].get_value(t)[r * static_cast<size_t>(mio::osir::InfectionState::Count) +
                                                static_cast<size_t>(mio::osir::InfectionState::Recovered)];
            }
        }

        std::array<int, static_cast<size_t>(mio::osir::InfectionState::Count) * NumRegions>
            indices; // Vector holding the
        std::function<void(int, int)> fill_moments =
            [&](int pos, int currentSum) { // pos: current position in indices, currentSum: sum of indices so far
                if (pos == int(indices.size()) &&
                    std::accumulate(indices.begin(), indices.end(), 0) <=
                        Order) { // Position is at last index i.e. all indiced for the moment are filled
                    moments[indices] = calculate_moment<NumRegions>(result, indices);
                    return;
                }

                int maxAllowedHere = std::min(Order, Order - currentSum); //maximum allowed value for current index
                for (int v = 0; v <= maxAllowedHere; ++v) { // Iterate over all values allowed for the current index
                    indices[pos] = v;
                    int newSum   = currentSum + v; // Increase sum by current index
                    fill_moments(
                        pos + 1,
                        newSum); // This triggers the next index to take all possible values given the value of the current index
                }
            };

        fill_moments(0, 0); // Start with first index and sum 0
        // Get only moments up to the given order
        auto moment_values   = moments.moments_up_to_order(Order);
        Eigen::VectorXd data = Eigen::VectorXd::Map(moment_values.data(), moment_values.size());
        moment_ts.add_time_point(sim_results[0].get_time(t), data);
    }
    return std::make_pair(moment_ts, moments.names_up_to_order(Order));
}

#endif // SMM_HELPER_H
