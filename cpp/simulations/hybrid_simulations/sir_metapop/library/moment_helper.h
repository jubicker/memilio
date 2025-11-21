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
#ifndef MOMENT_HELPER_H
#define MOMENT_HELPER_H

#include "hybrid_simulations/sir_metapop/config/config.h"
#include "hybrid_simulations/sir_metapop/library/moments/model.h"
#include <cstddef>

namespace moment_helper
{

template <size_t NumRegions, size_t ClosureOrder>
smm_moments::Model<NumRegions, ClosureOrder>
initialize_model(Eigen::Array<double, Eigen::Dynamic, 1>& expected_values_init,
                 Eigen::Array<double, Eigen::Dynamic, 1>& moments_init, const Config::Config& config)
{
    smm_moments::Model<NumRegions, ClosureOrder> model;
    assert(expected_values_init.rows() == NumRegions * static_cast<size_t>(mio::osir::InfectionState::Count) &&
           "Initial expected values do not have correct size");
    assert(moments_init.rows() == model.moments.moments().rows() && "Initial moments do not have correct size");
    // Set transmission rates
    for (auto& rate : config.transition_rates) {
        model.parameters.template get<smm_moments::TransmissionRate>()[{rate.status, rate.from, rate.to}] = rate.factor;
    }
    // Set recovery rate
    model.parameters.template get<smm_moments::RecoveryRate>() = config.gamma;
    for (size_t r = 0; r < NumRegions; ++r) {
        // Set initial expected values
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
    // Set initial moments
    model.moments.moments() = moments_init;
    return model;
}

} // namespace moment_helper
#endif // MOMENT_HELPER_H
