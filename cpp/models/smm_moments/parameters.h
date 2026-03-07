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

#ifndef MOMENTS_PARAMETERS_H
#define MOMENTS_PARAMETERS_H

#include "memilio/config.h"
#include "memilio/geography/regions.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/utils/parameter_set.h"
#include "ode_sir/infection_state.h"

namespace mio
{

namespace smm_moments
{

/**
 * @brief Rate from S to I. Is region-dependent.
 */
struct TransmissionRate {
    using Type = mio::CustomIndexArray<ScalarType, mio::regions::Region>;
    static Type get_default(mio::regions::Region size)
    {
        return Type(size, 1.0);
    }
    static std::string name()
    {
        return "TransmissionRate";
    }
};

/**
 * @brief Rate from I to R.
 */
struct RecoveryRate {
    using Type = mio::CustomIndexArray<ScalarType, mio::regions::Region>;
    static Type get_default(mio::regions::Region size)
    {
        return Type(size, 1.0);
    }
    static std::string name()
    {
        return "RecoveryRate";
    }
};

/**
 * @brief Spatial transition rate from region k to region l. Is dependent on infection state.
 */
struct TransitionRate {
    using Type =
        mio::CustomIndexArray<ScalarType, mio::osir::InfectionState, mio::regions::Region, mio::regions::Region>;

    static Type get_default(mio::regions::Region size)
    {
        return Type({mio::osir::InfectionState::Count, size, size}, 0.0);
    }
    static std::string name()
    {
        return "TransitionRate";
    }
};

using ParametersBase = mio::ParameterSet<TransmissionRate, RecoveryRate, TransitionRate>;

} // namespace smm_moments

} // namespace mio

#endif // MOMENTS_PARAMETERS_H
