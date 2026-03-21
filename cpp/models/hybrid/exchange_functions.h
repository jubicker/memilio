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

#ifndef MIO_HYBRID_EXCHANGE_FUNCTIONS_H
#define MIO_HYBRID_EXCHANGE_FUNCTIONS_H

#include "smm/simulation_set.h"
#include "ode_sir/infection_state.h"
#include "smm_moments/simulation.h"

namespace mio
{
namespace hybrid
{

template <class ModelFrom, class ModelTo>
void exchange_agents(ModelFrom&, ModelTo&, size_t, size_t) = delete;

template <>
void exchange_agents(mio::smm::SimulationSet<2, mio::osir::InfectionState, 3>& model_from,
                     mio::smm::SimulationSet<2, mio::osir::InfectionState, 3>& model_to, size_t /*region_from*/,
                     size_t region_to);

template <>
void exchange_agents(mio::smm_moments::Simulation<2, 3>& model_from, mio::smm_moments::Simulation<2, 3>& model_to,
                     size_t /*region_from*/, size_t region_to);

template <>
void exchange_agents(mio::smm::SimulationSet<2, mio::osir::InfectionState, 3>& model_from,
                     mio::smm_moments::Simulation<2, 3>& model_to, size_t /*region_from*/, size_t region_to);

template <>
void exchange_agents(mio::smm_moments::Simulation<2, 3>& model_from,
                     mio::smm::SimulationSet<2, mio::osir::InfectionState, 3>& model_to, size_t /*region_from*/,
                     size_t region_to);

} //namespace hybrid

} //namespace mio

#endif //MIO_HYBRID_EXCHANGE_FUNCTIONS_H
