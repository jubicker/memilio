/* 
* Copyright (C) 2020-2025 German Aerospace Center (DLR-SC)
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

#ifndef MOMENTS_MODEL_H
#define MOMENTS_MODEL_H

#include "memilio/config.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "simulations/hybrid_simulations/sir_metapop/library/moment_array.h"
#include "smm_moments/parameters.h"
#include "smm_moments/closure_functions.h"
#include <boost/math/special_functions/math_fwd.hpp>
#include <cstddef>

#include <boost/math/special_functions/binomial.hpp>

template double boost::math::binomial_coefficient<double>(unsigned, unsigned);
namespace mio
{

namespace smm_moments
{

/**
 * @brief Moment equation model for SIR-SMM.
 * @tparam NumRegions Number of regions.
 * @tparam ClosureOrder Order used for zero cumulant closure.
 */
template <size_t NumRegions, size_t ClosureOrder>
class Model
{
public:
    using Region         = mio::regions::Region;
    using InfectionState = mio::osir::InfectionState;

    using ClosureFunctionType = ScalarType (*)(
        size_t closure_order, std::array<int, static_cast<size_t>(osir::InfectionState::Count) * NumRegions> index,
        Eigen::Ref<const Eigen::VectorX<ScalarType>> y);

    Model(ClosureFunctionType closure_func = &truncation_closure<NumRegions>)
        : parameters(NumRegions)
        , populations({static_cast<Region>(NumRegions), InfectionState::Count}, 0.0)
        , m_closure_function(closure_func)
    {
    }

    /**
     * @brief Evaluates right-hand-side of ODEs describing expected values and moments up to ClosureOrder.
     * @param[in] y Current state of the expected values (first entries) and the moments (following entries).
     * @param[in] dydt Reference to the calculated output.
     */
    void get_derivatives(Eigen::Ref<const Eigen::VectorX<ScalarType>> y, ScalarType /*t*/,
                         Eigen::Ref<Eigen::VectorX<ScalarType>> dydt) const
    {
        std::array<int, static_cast<size_t>(InfectionState::Count) * NumRegions> indices;
        for (size_t l = 0; l < NumRegions; ++l) {
            indices.fill(0);
            // Indices for S, I, R in region l
            size_t Sl = this->populations.get_flat_index({Region(l), InfectionState::Susceptible});
            size_t Il = this->populations.get_flat_index({Region(l), InfectionState::Infected});
            size_t Rl = this->populations.get_flat_index({Region(l), InfectionState::Recovered});
            // Multiindex for M_1SlIl
            indices[l * static_cast<size_t>(InfectionState::Count) + static_cast<size_t>(InfectionState::Susceptible)] =
                1;
            indices[l * static_cast<size_t>(InfectionState::Count) + static_cast<size_t>(InfectionState::Infected)] = 1;
            assert(ClosureOrder >= 2);
            double M_1SlIl =
                ClosureOrder > 2 ? y[moments.flatten_index(indices) + populations.get_num_compartments()] : 0;
            dydt[Sl] = -parameters.template get<TransmissionRate>()[Region(l)] * (y[Sl] * y[Il] + M_1SlIl);
            dydt[Il] = parameters.template get<TransmissionRate>()[Region(l)] * (y[Sl] * y[Il] + M_1SlIl) -
                       parameters.template get<RecoveryRate>() * y[Il];
            dydt[Rl] = parameters.template get<RecoveryRate>() * y[Il];
            for (size_t k = 0; k < NumRegions; ++k) {
                if (k == l) {
                    continue;
                }
                // Indices for S, I, R in region k
                size_t Sk = this->populations.get_flat_index({Region(k), InfectionState::Susceptible});
                size_t Ik = this->populations.get_flat_index({Region(k), InfectionState::Infected});
                size_t Rk = this->populations.get_flat_index({Region(k), InfectionState::Recovered});
                dydt[Sl] +=
                    parameters.template get<TransitionRate>()[{InfectionState::Susceptible, Region(k), Region(l)}] *
                        y[Sk] -
                    parameters.template get<TransitionRate>()[{InfectionState::Susceptible, Region(l), Region(k)}] *
                        y[Sl];
                dydt[Il] +=
                    parameters.template get<TransitionRate>()[{InfectionState::Infected, Region(k), Region(l)}] *
                        y[Ik] -
                    parameters.template get<TransitionRate>()[{InfectionState::Infected, Region(l), Region(k)}] * y[Il];
                dydt[Rl] +=
                    parameters.template get<TransitionRate>()[{InfectionState::Recovered, Region(k), Region(l)}] *
                        y[Rk] -
                    parameters.template get<TransitionRate>()[{InfectionState::Recovered, Region(l), Region(k)}] *
                        y[Rl];
            }
        }

        if (ClosureOrder > 2) {
            indices.fill(0);
            // recursive function that claculates derivatives for all moments up to order ClosureOrder
            std::function<void(int, int)> fill_moments =
                [&](int pos, int currentSum) { // pos: current position in indices, currentSum: sum of indices so far
                    if (pos == int(indices.size()) &&
                        std::accumulate(indices.begin(), indices.end(), 0) <
                            int(ClosureOrder)) { // Position is at last index i.e. all indiced for the moment are filled
                        if (std::accumulate(indices.begin(), indices.end(), 0) == 0 ||
                            std::accumulate(indices.begin(), indices.end(), 0) == 1) {
                            return; // Moment of order 0 is always 1 and moments of order 1 are always 0
                        }
                        size_t flat_index = moments.flatten_index(indices) +
                                            populations.get_num_compartments(); // flat index of current moment
                        dydt[flat_index] = 0.;
                        for (size_t l = 0; l < NumRegions; ++l) {
                            // Indices for S, I, R in region l
                            size_t Sl    = this->populations.get_flat_index({Region(l), InfectionState::Susceptible});
                            size_t Il    = this->populations.get_flat_index({Region(l), InfectionState::Infected});
                            size_t Rl    = this->populations.get_flat_index({Region(l), InfectionState::Recovered});
                            size_t i_S_l = indices[l * static_cast<size_t>(InfectionState::Count) +
                                                   static_cast<size_t>(InfectionState::Susceptible)];
                            size_t i_I_l = indices[l * static_cast<size_t>(InfectionState::Count) +
                                                   static_cast<size_t>(InfectionState::Infected)];
                            size_t i_R_l = indices[l * static_cast<size_t>(InfectionState::Count) +
                                                   static_cast<size_t>(InfectionState::Recovered)];
                            for (size_t h_S_l = 0; h_S_l <= i_S_l; ++h_S_l) {
                                for (size_t h_I_l = 0; h_I_l <= i_I_l; ++h_I_l) {
                                    if (h_S_l + h_I_l == i_S_l + i_I_l) {
                                        continue;
                                    }
                                    std::array<int, static_cast<size_t>(InfectionState::Count)* NumRegions>
                                        h_S_l_h_I_l_index                                               = indices;
                                    h_S_l_h_I_l_index[l * static_cast<size_t>(InfectionState::Count) +
                                                      static_cast<size_t>(InfectionState::Susceptible)] = h_S_l;
                                    h_S_l_h_I_l_index[l * static_cast<size_t>(InfectionState::Count) +
                                                      static_cast<size_t>(InfectionState::Infected)]    = h_I_l;
                                    int current_order =
                                        std::accumulate(h_S_l_h_I_l_index.begin(), h_S_l_h_I_l_index.end(), 0);
                                    double M_h_S_l_h_I_l = int(ClosureOrder) > current_order
                                                               ? y[moments.flatten_index(h_S_l_h_I_l_index) +
                                                                   populations.get_num_compartments()]
                                                               : m_closure_function(ClosureOrder, h_S_l_h_I_l_index, y);
                                    if (current_order == 0) {
                                        M_h_S_l_h_I_l = 1.;
                                    }
                                    dydt[flat_index] += parameters.template get<TransmissionRate>()[Region(l)] * y[Sl] *
                                                        y[Il] *
                                                        boost::math::binomial_coefficient<double>(i_S_l, h_S_l) *
                                                        boost::math::binomial_coefficient<double>(i_I_l, h_I_l) *
                                                        std::pow(-1, i_S_l - h_S_l) * M_h_S_l_h_I_l;
                                } // h_I_l
                            } // h_S_l

                            for (size_t h_I_l = 0; h_I_l <= i_I_l; ++h_I_l) {
                                for (size_t h_R_l = 0; h_R_l <= i_R_l; ++h_R_l) {
                                    if (h_I_l + h_R_l == i_I_l + i_R_l) {
                                        continue;
                                    }
                                    std::array<int, static_cast<size_t>(InfectionState::Count)* NumRegions>
                                        h_I_l_h_R_l_index                                             = indices;
                                    h_I_l_h_R_l_index[l * static_cast<size_t>(InfectionState::Count) +
                                                      static_cast<size_t>(InfectionState::Infected)]  = h_I_l;
                                    h_I_l_h_R_l_index[l * static_cast<size_t>(InfectionState::Count) +
                                                      static_cast<size_t>(InfectionState::Recovered)] = h_R_l;
                                    int current_order =
                                        std::accumulate(h_I_l_h_R_l_index.begin(), h_I_l_h_R_l_index.end(), 0);
                                    double M_h_I_l_h_R_l = int(ClosureOrder) > current_order
                                                               ? y[moments.flatten_index(h_I_l_h_R_l_index) +
                                                                   populations.get_num_compartments()]
                                                               : m_closure_function(ClosureOrder, h_I_l_h_R_l_index, y);
                                    if (current_order == 0) {
                                        M_h_I_l_h_R_l = 1.;
                                    }
                                    dydt[flat_index] += parameters.template get<RecoveryRate>() * y[Il] *
                                                        boost::math::binomial_coefficient<double>(i_I_l, h_I_l) *
                                                        boost::math::binomial_coefficient<double>(i_R_l, h_R_l) *
                                                        std::pow(-1, i_I_l - h_I_l) * M_h_I_l_h_R_l;
                                } // h_R_l
                            } // h_I_l

                            for (size_t k = l + 1; k < NumRegions; ++k) {
                                size_t Sk = this->populations.get_flat_index({Region(k), InfectionState::Susceptible});
                                size_t Ik = this->populations.get_flat_index({Region(k), InfectionState::Infected});
                                size_t Rk = this->populations.get_flat_index({Region(k), InfectionState::Recovered});
                                size_t i_S_k = indices[k * static_cast<size_t>(InfectionState::Count) +
                                                       static_cast<size_t>(InfectionState::Susceptible)];
                                size_t i_I_k = indices[k * static_cast<size_t>(InfectionState::Count) +
                                                       static_cast<size_t>(InfectionState::Infected)];
                                size_t i_R_k = indices[k * static_cast<size_t>(InfectionState::Count) +
                                                       static_cast<size_t>(InfectionState::Recovered)];
                                for (size_t h_S_l = 0; h_S_l <= i_S_l; ++h_S_l) {
                                    for (size_t h_S_k = 0; h_S_k <= i_S_k; ++h_S_k) {
                                        if (h_S_l + h_S_k == i_S_l + i_S_k) {
                                            continue;
                                        }
                                        std::array<int, static_cast<size_t>(InfectionState::Count)* NumRegions>
                                            h_S_l_h_S_k_index                                               = indices;
                                        h_S_l_h_S_k_index[l * static_cast<size_t>(InfectionState::Count) +
                                                          static_cast<size_t>(InfectionState::Susceptible)] = h_S_l;
                                        h_S_l_h_S_k_index[k * static_cast<size_t>(InfectionState::Count) +
                                                          static_cast<size_t>(InfectionState::Susceptible)] = h_S_k;
                                        int current_order =
                                            std::accumulate(h_S_l_h_S_k_index.begin(), h_S_l_h_S_k_index.end(), 0);
                                        double M_h_S_l_h_S_k =
                                            int(ClosureOrder) > current_order
                                                ? y[moments.flatten_index(h_S_l_h_S_k_index) +
                                                    populations.get_num_compartments()]
                                                : m_closure_function(ClosureOrder, h_S_l_h_S_k_index, y);
                                        if (current_order == 0) {
                                            M_h_S_l_h_S_k = 1.;
                                        }
                                        dydt[flat_index] +=
                                            boost::math::binomial_coefficient<double>(i_S_l, h_S_l) *
                                            boost::math::binomial_coefficient<double>(i_S_k, h_S_k) * M_h_S_l_h_S_k *
                                            (std::pow(-1, i_S_l - h_S_l) *
                                                 parameters.template get<TransitionRate>()[{InfectionState::Susceptible,
                                                                                            Region(l), Region(k)}] *
                                                 y[Sl] +
                                             std::pow(-1, i_S_k - h_S_k) *
                                                 parameters.template get<TransitionRate>()[{InfectionState::Susceptible,
                                                                                            Region(k), Region(l)}] *
                                                 y[Sk]);
                                    } // h_S_k
                                } // h_S_l

                                for (size_t h_I_l = 0; h_I_l <= i_I_l; ++h_I_l) {
                                    for (size_t h_I_k = 0; h_I_k <= i_I_k; ++h_I_k) {
                                        if (h_I_l + h_I_k == i_I_l + i_I_k) {
                                            continue;
                                        }
                                        std::array<int, static_cast<size_t>(InfectionState::Count)* NumRegions>
                                            h_I_l_h_I_k_index                                            = indices;
                                        h_I_l_h_I_k_index[l * static_cast<size_t>(InfectionState::Count) +
                                                          static_cast<size_t>(InfectionState::Infected)] = h_I_l;
                                        h_I_l_h_I_k_index[k * static_cast<size_t>(InfectionState::Count) +
                                                          static_cast<size_t>(InfectionState::Infected)] = h_I_k;
                                        int current_order =
                                            std::accumulate(h_I_l_h_I_k_index.begin(), h_I_l_h_I_k_index.end(), 0);
                                        double M_h_I_l_h_I_k =
                                            int(ClosureOrder) > current_order
                                                ? y[moments.flatten_index(h_I_l_h_I_k_index) +
                                                    populations.get_num_compartments()]
                                                : m_closure_function(ClosureOrder, h_I_l_h_I_k_index, y);
                                        if (current_order == 0) {
                                            M_h_I_l_h_I_k = 1.;
                                        }
                                        dydt[flat_index] +=
                                            boost::math::binomial_coefficient<double>(i_I_l, h_I_l) *
                                            boost::math::binomial_coefficient<double>(i_I_k, h_I_k) * M_h_I_l_h_I_k *
                                            (std::pow(-1, i_I_l - h_I_l) *
                                                 parameters.template get<TransitionRate>()[{InfectionState::Infected,
                                                                                            Region(l), Region(k)}] *
                                                 y[Il] +
                                             std::pow(-1, i_I_k - h_I_k) *
                                                 parameters.template get<TransitionRate>()[{InfectionState::Infected,
                                                                                            Region(k), Region(l)}] *
                                                 y[Ik]);
                                    } // h_I_k
                                } // h_I_l

                                for (size_t h_R_l = 0; h_R_l <= i_R_l; ++h_R_l) {
                                    for (size_t h_R_k = 0; h_R_k <= i_R_k; ++h_R_k) {
                                        if (h_R_l + h_R_k == i_R_l + i_R_k) {
                                            continue;
                                        }
                                        std::array<int, static_cast<size_t>(InfectionState::Count)* NumRegions>
                                            h_R_l_h_R_k_index                                             = indices;
                                        h_R_l_h_R_k_index[l * static_cast<size_t>(InfectionState::Count) +
                                                          static_cast<size_t>(InfectionState::Recovered)] = h_R_l;
                                        h_R_l_h_R_k_index[k * static_cast<size_t>(InfectionState::Count) +
                                                          static_cast<size_t>(InfectionState::Recovered)] = h_R_k;
                                        int current_order =
                                            std::accumulate(h_R_l_h_R_k_index.begin(), h_R_l_h_R_k_index.end(), 0);
                                        double M_h_R_l_h_R_k =
                                            int(ClosureOrder) > current_order
                                                ? y[moments.flatten_index(h_R_l_h_R_k_index) +
                                                    populations.get_num_compartments()]
                                                : m_closure_function(ClosureOrder, h_R_l_h_R_k_index, y);
                                        if (current_order == 0) {
                                            M_h_R_l_h_R_k = 1.;
                                        }
                                        dydt[flat_index] +=
                                            boost::math::binomial_coefficient<double>(i_R_l, h_R_l) *
                                            boost::math::binomial_coefficient<double>(i_R_k, h_R_k) * M_h_R_l_h_R_k *
                                            (std::pow(-1, i_R_l - h_R_l) *
                                                 parameters.template get<TransitionRate>()[{InfectionState::Recovered,
                                                                                            Region(l), Region(k)}] *
                                                 y[Rl] +
                                             std::pow(-1, i_R_k - h_R_k) *
                                                 parameters.template get<TransitionRate>()[{InfectionState::Recovered,
                                                                                            Region(k), Region(l)}] *
                                                 y[Rk]);
                                    } // h_R_k
                                } // h_R_l
                            } // k
                            // First derivatives
                            for (size_t h_S_l = 0; h_S_l <= i_S_l; ++h_S_l) {
                                for (size_t h_I_l = 0; h_I_l <= i_I_l; ++h_I_l) {
                                    if (h_S_l + h_I_l == i_S_l + i_I_l) {
                                        continue;
                                    }
                                    std::array<int, static_cast<size_t>(InfectionState::Count)* NumRegions>
                                        h_S_l_p1_h_I_l_index                                               = indices;
                                    h_S_l_p1_h_I_l_index[l * static_cast<size_t>(InfectionState::Count) +
                                                         static_cast<size_t>(InfectionState::Susceptible)] = h_S_l + 1;
                                    h_S_l_p1_h_I_l_index[l * static_cast<size_t>(InfectionState::Count) +
                                                         static_cast<size_t>(InfectionState::Infected)]    = h_I_l;
                                    int current_order1 =
                                        std::accumulate(h_S_l_p1_h_I_l_index.begin(), h_S_l_p1_h_I_l_index.end(), 0);
                                    double M_h_S_l_p1_h_I_l =
                                        int(ClosureOrder) > current_order1
                                            ? y[moments.flatten_index(h_S_l_p1_h_I_l_index) +
                                                populations.get_num_compartments()]
                                            : m_closure_function(ClosureOrder, h_S_l_p1_h_I_l_index, y);
                                    if (current_order1 == 0) {
                                        M_h_S_l_p1_h_I_l = 1.;
                                    }
                                    std::array<int, static_cast<size_t>(InfectionState::Count)* NumRegions>
                                        h_S_l_h_I_l_p1_index                                               = indices;
                                    h_S_l_h_I_l_p1_index[l * static_cast<size_t>(InfectionState::Count) +
                                                         static_cast<size_t>(InfectionState::Susceptible)] = h_S_l;
                                    h_S_l_h_I_l_p1_index[l * static_cast<size_t>(InfectionState::Count) +
                                                         static_cast<size_t>(InfectionState::Infected)]    = h_I_l + 1;
                                    int current_order2 =
                                        std::accumulate(h_S_l_h_I_l_p1_index.begin(), h_S_l_h_I_l_p1_index.end(), 0);
                                    double M_h_S_l_h_I_l_p1 =
                                        int(ClosureOrder) > current_order2
                                            ? y[moments.flatten_index(h_S_l_h_I_l_p1_index) +
                                                populations.get_num_compartments()]
                                            : m_closure_function(ClosureOrder, h_S_l_h_I_l_p1_index, y);
                                    if (current_order2 == 0) {
                                        M_h_S_l_h_I_l_p1 = 1.;
                                    }
                                    dydt[flat_index] += boost::math::binomial_coefficient<double>(i_S_l, h_S_l) *
                                                        boost::math::binomial_coefficient<double>(i_I_l, h_I_l) *
                                                        std::pow(-1, i_S_l - h_S_l) *
                                                        parameters.template get<TransmissionRate>()[Region(l)] *
                                                        (y[Il] * M_h_S_l_p1_h_I_l + y[Sl] * M_h_S_l_h_I_l_p1);

                                } // h_I_l
                            } // h_S_l

                            for (size_t h_I_l = 0; h_I_l <= i_I_l; ++h_I_l) {
                                for (size_t h_R_l = 0; h_R_l <= i_R_l; ++h_R_l) {
                                    if (h_I_l + h_R_l == i_I_l + i_R_l) {
                                        continue;
                                    }
                                    std::array<int, static_cast<size_t>(InfectionState::Count)* NumRegions>
                                        h_I_l_p1_h_R_l_index                                             = indices;
                                    h_I_l_p1_h_R_l_index[l * static_cast<size_t>(InfectionState::Count) +
                                                         static_cast<size_t>(InfectionState::Infected)]  = h_I_l + 1;
                                    h_I_l_p1_h_R_l_index[l * static_cast<size_t>(InfectionState::Count) +
                                                         static_cast<size_t>(InfectionState::Recovered)] = h_R_l;
                                    int current_order =
                                        std::accumulate(h_I_l_p1_h_R_l_index.begin(), h_I_l_p1_h_R_l_index.end(), 0);
                                    double M_h_I_l_p1_h_R_l =
                                        int(ClosureOrder) > current_order
                                            ? y[moments.flatten_index(h_I_l_p1_h_R_l_index) +
                                                populations.get_num_compartments()]
                                            : m_closure_function(ClosureOrder, h_I_l_p1_h_R_l_index, y);
                                    if (current_order == 0) {
                                        M_h_I_l_p1_h_R_l = 1.;
                                    }
                                    dydt[flat_index] += parameters.template get<RecoveryRate>() *
                                                        boost::math::binomial_coefficient<double>(i_I_l, h_I_l) *
                                                        boost::math::binomial_coefficient<double>(i_R_l, h_R_l) *
                                                        std::pow(-1, i_I_l - h_I_l) * M_h_I_l_p1_h_R_l;
                                } // h_R_l
                            } // h_I_l

                            for (size_t k = l + 1; k < NumRegions; ++k) {
                                size_t i_S_k = indices[k * static_cast<size_t>(InfectionState::Count) +
                                                       static_cast<size_t>(InfectionState::Susceptible)];
                                size_t i_I_k = indices[k * static_cast<size_t>(InfectionState::Count) +
                                                       static_cast<size_t>(InfectionState::Infected)];
                                size_t i_R_k = indices[k * static_cast<size_t>(InfectionState::Count) +
                                                       static_cast<size_t>(InfectionState::Recovered)];
                                for (size_t h_S_l = 0; h_S_l <= i_S_l; ++h_S_l) {
                                    for (size_t h_S_k = 0; h_S_k <= i_S_k; ++h_S_k) {
                                        if (h_S_l + h_S_k == i_S_l + i_S_k) {
                                            continue;
                                        }
                                        std::array<int, static_cast<size_t>(InfectionState::Count)* NumRegions>
                                            h_S_l_p1_h_S_k_index = indices;
                                        h_S_l_p1_h_S_k_index[l * static_cast<size_t>(InfectionState::Count) +
                                                             static_cast<size_t>(InfectionState::Susceptible)] =
                                            h_S_l + 1;
                                        h_S_l_p1_h_S_k_index[k * static_cast<size_t>(InfectionState::Count) +
                                                             static_cast<size_t>(InfectionState::Susceptible)] = h_S_k;
                                        int current_order1 = std::accumulate(h_S_l_p1_h_S_k_index.begin(),
                                                                             h_S_l_p1_h_S_k_index.end(), 0);
                                        double M_h_S_l_p1_h_S_k =
                                            int(ClosureOrder) > current_order1
                                                ? y[moments.flatten_index(h_S_l_p1_h_S_k_index) +
                                                    populations.get_num_compartments()]
                                                : m_closure_function(ClosureOrder, h_S_l_p1_h_S_k_index, y);
                                        if (current_order1 == 0) {
                                            M_h_S_l_p1_h_S_k = 1.;
                                        }
                                        std::array<int, static_cast<size_t>(InfectionState::Count)* NumRegions>
                                            h_S_l_h_S_k_p1_index = indices;
                                        h_S_l_h_S_k_p1_index[l * static_cast<size_t>(InfectionState::Count) +
                                                             static_cast<size_t>(InfectionState::Susceptible)] = h_S_l;
                                        h_S_l_h_S_k_p1_index[k * static_cast<size_t>(InfectionState::Count) +
                                                             static_cast<size_t>(InfectionState::Susceptible)] =
                                            h_S_k + 1;
                                        int current_order2 = std::accumulate(h_S_l_h_S_k_p1_index.begin(),
                                                                             h_S_l_h_S_k_p1_index.end(), 0);
                                        double M_h_S_l_h_S_k_p1 =
                                            int(ClosureOrder) > current_order2
                                                ? y[moments.flatten_index(h_S_l_h_S_k_p1_index) +
                                                    populations.get_num_compartments()]
                                                : m_closure_function(ClosureOrder, h_S_l_h_S_k_p1_index, y);
                                        if (current_order2 == 0) {
                                            M_h_S_l_h_S_k_p1 = 1.;
                                        }
                                        dydt[flat_index] +=
                                            boost::math::binomial_coefficient<double>(i_S_l, h_S_l) *
                                            boost::math::binomial_coefficient<double>(i_S_k, h_S_k) *
                                            (std::pow(-1, i_S_l - h_S_l) *
                                                 parameters.template get<TransitionRate>()[{InfectionState::Susceptible,
                                                                                            Region(l), Region(k)}] *
                                                 M_h_S_l_p1_h_S_k +
                                             std::pow(-1, i_S_k - h_S_k) *
                                                 parameters.template get<TransitionRate>()[{InfectionState::Susceptible,
                                                                                            Region(k), Region(l)}] *
                                                 M_h_S_l_h_S_k_p1);
                                    } // h_S_k
                                } // h_S_l

                                for (size_t h_I_l = 0; h_I_l <= i_I_l; ++h_I_l) {
                                    for (size_t h_I_k = 0; h_I_k <= i_I_k; ++h_I_k) {
                                        if (h_I_l + h_I_k == i_I_l + i_I_k) {
                                            continue;
                                        }
                                        std::array<int, static_cast<size_t>(InfectionState::Count)* NumRegions>
                                            h_I_l_p1_h_I_k_index                                            = indices;
                                        h_I_l_p1_h_I_k_index[l * static_cast<size_t>(InfectionState::Count) +
                                                             static_cast<size_t>(InfectionState::Infected)] = h_I_l + 1;
                                        h_I_l_p1_h_I_k_index[k * static_cast<size_t>(InfectionState::Count) +
                                                             static_cast<size_t>(InfectionState::Infected)] = h_I_k;
                                        int current_order1 = std::accumulate(h_I_l_p1_h_I_k_index.begin(),
                                                                             h_I_l_p1_h_I_k_index.end(), 0);
                                        double M_h_I_l_p1_h_I_k =
                                            int(ClosureOrder) > current_order1
                                                ? y[moments.flatten_index(h_I_l_p1_h_I_k_index) +
                                                    populations.get_num_compartments()]
                                                : m_closure_function(ClosureOrder, h_I_l_p1_h_I_k_index, y);
                                        if (current_order1 == 0) {
                                            M_h_I_l_p1_h_I_k = 1.;
                                        }
                                        std::array<int, static_cast<size_t>(InfectionState::Count)* NumRegions>
                                            h_I_l_h_I_k_p1_index                                            = indices;
                                        h_I_l_h_I_k_p1_index[l * static_cast<size_t>(InfectionState::Count) +
                                                             static_cast<size_t>(InfectionState::Infected)] = h_I_l;
                                        h_I_l_h_I_k_p1_index[k * static_cast<size_t>(InfectionState::Count) +
                                                             static_cast<size_t>(InfectionState::Infected)] = h_I_k + 1;
                                        int current_order2 = std::accumulate(h_I_l_h_I_k_p1_index.begin(),
                                                                             h_I_l_h_I_k_p1_index.end(), 0);
                                        double M_h_I_l_h_I_k_p1 =
                                            int(ClosureOrder) > current_order2
                                                ? y[moments.flatten_index(h_I_l_h_I_k_p1_index) +
                                                    populations.get_num_compartments()]
                                                : m_closure_function(ClosureOrder, h_I_l_h_I_k_p1_index, y);
                                        if (current_order2 == 0) {
                                            M_h_I_l_h_I_k_p1 = 1.;
                                        }
                                        dydt[flat_index] += boost::math::binomial_coefficient<double>(i_I_l, h_I_l) *
                                                            boost::math::binomial_coefficient<double>(i_I_k, h_I_k) *
                                                            (std::pow(-1, i_I_l - h_I_l) *
                                                                 parameters.template get<TransitionRate>()[{
                                                                     InfectionState::Infected, Region(l), Region(k)}] *
                                                                 M_h_I_l_p1_h_I_k +
                                                             std::pow(-1, i_I_k - h_I_k) *
                                                                 parameters.template get<TransitionRate>()[{
                                                                     InfectionState::Infected, Region(k), Region(l)}] *
                                                                 M_h_I_l_h_I_k_p1);
                                    } // h_I_k
                                } // h_I_l

                                for (size_t h_R_l = 0; h_R_l <= i_R_l; ++h_R_l) {
                                    for (size_t h_R_k = 0; h_R_k <= i_R_k; ++h_R_k) {
                                        if (h_R_l + h_R_k == i_R_l + i_R_k) {
                                            continue;
                                        }
                                        std::array<int, static_cast<size_t>(InfectionState::Count)* NumRegions>
                                            h_R_l_p1_h_R_k_index = indices;
                                        h_R_l_p1_h_R_k_index[l * static_cast<size_t>(InfectionState::Count) +
                                                             static_cast<size_t>(InfectionState::Recovered)] =
                                            h_R_l + 1;
                                        h_R_l_p1_h_R_k_index[k * static_cast<size_t>(InfectionState::Count) +
                                                             static_cast<size_t>(InfectionState::Recovered)] = h_R_k;
                                        int current_order1 = std::accumulate(h_R_l_p1_h_R_k_index.begin(),
                                                                             h_R_l_p1_h_R_k_index.end(), 0);
                                        double M_h_R_l_p1_h_R_k =
                                            int(ClosureOrder) > current_order1
                                                ? y[moments.flatten_index(h_R_l_p1_h_R_k_index) +
                                                    populations.get_num_compartments()]
                                                : m_closure_function(ClosureOrder, h_R_l_p1_h_R_k_index, y);
                                        if (current_order1 == 0) {
                                            M_h_R_l_p1_h_R_k = 1.;
                                        }
                                        std::array<int, static_cast<size_t>(InfectionState::Count)* NumRegions>
                                            h_R_l_h_R_k_p1_index                                             = indices;
                                        h_R_l_h_R_k_p1_index[l * static_cast<size_t>(InfectionState::Count) +
                                                             static_cast<size_t>(InfectionState::Recovered)] = h_R_l;
                                        h_R_l_h_R_k_p1_index[k * static_cast<size_t>(InfectionState::Count) +
                                                             static_cast<size_t>(InfectionState::Recovered)] =
                                            h_R_k + 1;
                                        int current_order2 = std::accumulate(h_R_l_h_R_k_p1_index.begin(),
                                                                             h_R_l_h_R_k_p1_index.end(), 0);
                                        double M_h_R_l_h_R_k_p1 =
                                            int(ClosureOrder) > current_order2
                                                ? y[moments.flatten_index(h_R_l_h_R_k_p1_index) +
                                                    populations.get_num_compartments()]
                                                : m_closure_function(ClosureOrder, h_R_l_h_R_k_p1_index, y);
                                        if (current_order2 == 0) {
                                            M_h_R_l_h_R_k_p1 = 1.;
                                        }
                                        dydt[flat_index] += boost::math::binomial_coefficient<double>(i_R_l, h_R_l) *
                                                            boost::math::binomial_coefficient<double>(i_R_k, h_R_k) *
                                                            (std::pow(-1, i_R_l - h_R_l) *
                                                                 parameters.template get<TransitionRate>()[{
                                                                     InfectionState::Recovered, Region(l), Region(k)}] *
                                                                 M_h_R_l_p1_h_R_k +
                                                             std::pow(-1, i_R_k - h_R_k) *
                                                                 parameters.template get<TransitionRate>()[{
                                                                     InfectionState::Recovered, Region(k), Region(l)}] *
                                                                 M_h_R_l_h_R_k_p1);
                                    } // h_R_k
                                } // h_R_l
                            } // k
                            // Second derivatives
                            for (size_t h_S_l = 0; h_S_l <= i_S_l; ++h_S_l) {
                                for (size_t h_I_l = 0; h_I_l <= i_I_l; ++h_I_l) {
                                    if (h_S_l + h_I_l == i_S_l + i_I_l) {
                                        continue;
                                    }
                                    std::array<int, static_cast<size_t>(InfectionState::Count)* NumRegions>
                                        h_S_l_p1_h_I_l_p1_index = indices;
                                    h_S_l_p1_h_I_l_p1_index[l * static_cast<size_t>(InfectionState::Count) +
                                                            static_cast<size_t>(InfectionState::Susceptible)] =
                                        h_S_l + 1;
                                    h_S_l_p1_h_I_l_p1_index[l * static_cast<size_t>(InfectionState::Count) +
                                                            static_cast<size_t>(InfectionState::Infected)] = h_I_l + 1;
                                    int current_order = std::accumulate(h_S_l_p1_h_I_l_p1_index.begin(),
                                                                        h_S_l_p1_h_I_l_p1_index.end(), 0);
                                    double M_h_S_l_p1_h_I_l_p1 =
                                        int(ClosureOrder) > current_order
                                            ? y[moments.flatten_index(h_S_l_p1_h_I_l_p1_index) +
                                                populations.get_num_compartments()]
                                            : m_closure_function(ClosureOrder, h_S_l_p1_h_I_l_p1_index, y);
                                    if (current_order == 0) {
                                        M_h_S_l_p1_h_I_l_p1 = 1.;
                                    }
                                    dydt[flat_index] += boost::math::binomial_coefficient<double>(i_S_l, h_S_l) *
                                                        boost::math::binomial_coefficient<double>(i_I_l, h_I_l) *
                                                        std::pow(-1, i_S_l - h_S_l) *
                                                        parameters.template get<TransmissionRate>()[Region(l)] *
                                                        M_h_S_l_p1_h_I_l_p1;

                                } // h_I_l
                            } // h_S_l
                            // Rest
                            std::array<int, static_cast<size_t>(InfectionState::Count)* NumRegions> i_S_l_m1_index =
                                indices;
                            i_S_l_m1_index[l * static_cast<size_t>(InfectionState::Count) +
                                           static_cast<size_t>(InfectionState::Susceptible)] = i_S_l - 1;
                            double M_i_S_l_m1                                                = 0;
                            if (i_S_l > 0) {
                                int current_order = std::accumulate(i_S_l_m1_index.begin(), i_S_l_m1_index.end(), 0);
                                M_i_S_l_m1 =
                                    int(ClosureOrder) > current_order
                                        ? y[moments.flatten_index(i_S_l_m1_index) + populations.get_num_compartments()]
                                        : m_closure_function(ClosureOrder, i_S_l_m1_index, y);
                                if (current_order == 0) {
                                    M_i_S_l_m1 = 1.;
                                }
                            }
                            dydt[flat_index] -= i_S_l * dydt[Sl] * M_i_S_l_m1;

                            std::array<int, static_cast<size_t>(InfectionState::Count)* NumRegions> i_I_l_m1_index =
                                indices;
                            i_I_l_m1_index[l * static_cast<size_t>(InfectionState::Count) +
                                           static_cast<size_t>(InfectionState::Infected)] = i_I_l - 1;
                            double M_i_I_l_m1                                             = 0;
                            if (i_I_l > 0) {
                                int current_order = std::accumulate(i_I_l_m1_index.begin(), i_I_l_m1_index.end(), 0);
                                M_i_I_l_m1 =
                                    int(ClosureOrder) > current_order
                                        ? y[moments.flatten_index(i_I_l_m1_index) + populations.get_num_compartments()]
                                        : m_closure_function(ClosureOrder, i_I_l_m1_index, y);
                                if (current_order == 0) {
                                    M_i_I_l_m1 = 1.;
                                }
                            }
                            dydt[flat_index] -= i_I_l * dydt[Il] * M_i_I_l_m1;

                            std::array<int, static_cast<size_t>(InfectionState::Count)* NumRegions> i_R_l_m1_index =
                                indices;
                            i_R_l_m1_index[l * static_cast<size_t>(InfectionState::Count) +
                                           static_cast<size_t>(InfectionState::Recovered)] = i_R_l - 1;
                            double M_i_R_l_m1                                              = 0;
                            if (i_R_l > 0) {
                                int current_order = std::accumulate(i_R_l_m1_index.begin(), i_R_l_m1_index.end(), 0);
                                M_i_R_l_m1 =
                                    int(ClosureOrder) > current_order
                                        ? y[moments.flatten_index(i_R_l_m1_index) + populations.get_num_compartments()]
                                        : m_closure_function(ClosureOrder, i_R_l_m1_index, y);
                                if (current_order == 0) {
                                    M_i_R_l_m1 = 1.;
                                }
                            }
                            dydt[flat_index] -= i_R_l * dydt[Rl] * M_i_R_l_m1;
                        } // l
                        return;
                    }

                    int maxAllowedHere = std::min(
                        ClosureOrder - 1, ClosureOrder - 1 - currentSum); //maximum allowed value for current index
                    for (int v = 0; v <= maxAllowedHere; ++v) { // Iterate over all values allowed for the current index
                        indices[pos] = v;
                        int newSum   = currentSum + v; // Increase sum by current index
                        fill_moments(
                            pos + 1,
                            newSum); // This triggers the next index to take all possible values given the value of the current index
                    }
                };
            fill_moments(0, 0); // Start with first index and sum 0
        }
    }

    /**
     * @brief Get the initial conditions for the ODE dydt = f(y, t).
     * See get_derivatives for more detail.
     * @return Current value of expected values and all moments as a flat vector.
     */
    Eigen::VectorX<ScalarType> get_initial_values() const
    {
        // Expected values are stored as Populations
        auto expected_values = populations.get_compartments();
        // Moment array (this contains also entries for moments > ClosureOrder as the MomentArray gets ClosureOrder as maximum value for every index)
        auto moment_vec = moments.moments();
        Eigen::VectorX<ScalarType> initial_values(expected_values.size() + moment_vec.size());
        initial_values.setZero();
        // Expected values are the first vector entries
        initial_values.head(expected_values.size()) = expected_values;
        // Moment values are copied to the vector afterwards
        for (auto i = 0; i < moment_vec.size(); ++i) {
            initial_values[expected_values.size() + i] = moment_vec[i];
        }
        std::array<int, static_cast<size_t>(InfectionState::Count) * NumRegions> zero_index;
        zero_index.fill(0);
        // M0 is 1
        initial_values[expected_values.size() + moments.flatten_index(zero_index)] = 1.0;
        return initial_values;
    }

    /**
     * @brief This function evaluates the right-hand-side f of the ODE dydt = f(y, t).
     * See get_derivatives.
     */
    void eval_right_hand_side(Eigen::Ref<const Eigen::VectorX<ScalarType>> /*pop*/,
                              Eigen::Ref<const Eigen::VectorX<ScalarType>> y, ScalarType t,
                              Eigen::Ref<Eigen::VectorX<ScalarType>> dydt) const
    {
        dydt.setZero();
        this->get_derivatives(y, t, dydt);
    }

    ParametersBase parameters{}; ///< Model's parameter set.
    MomentArray<static_cast<size_t>(InfectionState::Count), NumRegions, ClosureOrder>
        moments{}; ///< Array with initial moment values
    mio::Populations<ScalarType, Region, InfectionState> populations; ///< Array with initial values for expected values

private:
    ClosureFunctionType m_closure_function; ///< Moment closure approximation function
};

} // namespace smm_moments

} // namespace mio

#endif // MOMENTS_MODEL_H
