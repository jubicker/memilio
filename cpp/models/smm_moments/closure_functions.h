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

#ifndef MOMENTS_CLOSURE_FUNCTIONS_H
#define MOMENTS_CLOSURE_FUNCTIONS_H

#include "memilio/config.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "models/ode_sir/infection_state.h"
#include "memilio/math/eigen.h"
#include <array>
#include <cstddef>

namespace mio
{
namespace smm_moments
{

template <size_t NumRegions>
ScalarType truncation_closure(size_t closure_order,
                              std::array<int, static_cast<size_t>(osir::InfectionState::Count) * NumRegions> index,
                              Eigen::Ref<const Eigen::VectorX<ScalarType>> y)
{
    mio::unused(closure_order, index, y);
    return 0;
}

} // namespace smm_moments
} // namespace mio

#endif // MOMENTS_CLOSURE_FUNCTIONS_H
