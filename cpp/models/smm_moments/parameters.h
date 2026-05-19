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

struct ImmunityLossRate {
    using Type = mio::CustomIndexArray<ScalarType, mio::regions::Region>;
    static Type get_default(mio::regions::Region size)
    {
        return Type(size, 0.0);
    }
    static std::string name()
    {
        return "ImmunityLossRate";
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

/**
 * @brief The standard deviation of the seasonality in the current season, used for calculating the seasonality factor..
 */
struct SeasonalitySigma {
    using Type = std::vector<double>;
    static Type get_default()
    {
        return Type(0, 50.0);
    }
    const static std::string name()
    {
        return "SeasonalitySigma";
    }
};

/**
 * @brief The start day of the first season, used for calculating the seasonality factor.
 * The day is given in days from 1st January in the first season, e.g. 0 for 1st January, -185 for 1st July (in the year before 1st January of that season), etc.
 */
struct FirstSeasonStartDay {
    using Type = int;
    static Type get_default()
    {
        return Type(0);
    }
    const static std::string name()
    {
        return "FirstSeasonStartDay";
    }
};

/**
 * @brief The start day in the model.
 * The start day is given in days from 1st January of the current season e.g. 0 for 1st January, 180 for 1st July, etc..
 * The day is given in days from 1st January in the first season, e.g. 0 for 1st January, -185 for 1st July (in the year before 1st January of that season), etc.
 */
struct StartDay {
    using Type = int;
    static Type get_default()
    {
        return Type(0);
    }
    const static std::string name()
    {
        return "StartDay";
    }
};

/**
 * @brief The day of the year with the peak transmission rate in the current season, used for calculating the seasonality factor.
 * The day is given in days from 1st January of the current season, e.g. 0 for 1st January, 180 for 1st July, etc.
 */
struct SeasonalityPeak {
    using Type = std::vector<int>;
    static Type get_default()
    {
        return Type(0, 0.);
    }
    const static std::string name()
    {
        return "SeasonalityPeak";
    }
};

/**
 * @brief The strength of the seasonality in the current season, used for calculating the seasonality factor.
 * The value is given as a factor between 0 and 1, where 1 means no variation and 0 means full variation.
 */
struct SeasonalityRho {
    using Type = std::vector<double>;
    static Type get_default()
    {
        return Type(0, 1.0);
    }
    const static std::string name()
    {
        return "SeasonalityRho";
    }
};

using ParametersBase =
    mio::ParameterSet<TransmissionRate, RecoveryRate, ImmunityLossRate, TransitionRate, SeasonalitySigma,
                      FirstSeasonStartDay, StartDay, SeasonalityPeak, SeasonalityRho>;

} // namespace smm_moments

} // namespace mio

#endif // MOMENTS_PARAMETERS_H
