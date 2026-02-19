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

#include "memilio/config.h"
#include "models/ode_sir/infection_state.h"
#include "smm_moments/closure_functions.h"
#include "memilio/math/eigen.h"
#include <array>
#include <cstddef>

namespace mio
{
namespace smm_moments
{

ScalarType get_log_mu(ScalarType var, ScalarType mean)
{
    return std::log(mean * mean / std::sqrt(var + mean * mean));
}

ScalarType get_log_sigma(ScalarType var, ScalarType mean)
{
    return std::sqrt(std::log(var / (mean * mean) + 1));
}

ScalarType get_log_cov(ScalarType cov, ScalarType mean1, ScalarType mean2)
{
    return std::log(1 + cov / (mean1 * mean2));
}

} // namespace smm_moments
} // namespace mio
