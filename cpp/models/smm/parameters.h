/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: René Schmieding, Julia Bicker
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

#ifndef MIO_SMM_PARAMETERS_H
#define MIO_SMM_PARAMETERS_H

#include "memilio/config.h"
#include "memilio/geography/regions.h"
#include "memilio/utils/parameter_set.h"
#include "memilio/epidemiology/adoption_rate.h"
#include <cstddef>
#include <vector>

namespace mio
{
namespace smm
{

/**
 * @brief A vector of AdoptionRate%s, see mio::AdoptionRate
 * @tparam Status An infection state enum.
 */
template <typename FP, class Status>
struct AdoptionRates {
    using Type = std::vector<AdoptionRate<FP, Status>>;
    const static std::string name()
    {
        return "AdoptionRates";
    }
};

/**
 * @brief Struct defining a possible regional transition in a Model based on Poisson Processes.
 * @tparam Status An infection state enum.
 */
template <typename FP, class Status>
struct TransitionRate {
    Status status; // i
    mio::regions::Region from; // k
    mio::regions::Region to; // l
    FP factor; // lambda_i^{kl}
};

template <typename FP, class Status>
struct TransitionRates {
    using Type = std::vector<TransitionRate<FP, Status>>;
    const static std::string name()
    {
        return "TransitionRates";
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
    static Type get_default(size_t size)
    {
        return Type(size, 0.);
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
    static Type get_default(size_t size)
    {
        return Type(size, 1.0);
    }
    const static std::string name()
    {
        return "SeasonalityRho";
    }
};

/**
 * @brief The standard deviation of the seasonality in the current season, used for calculating the seasonality factor..
 */
struct SeasonalitySigma {
    using Type = std::vector<double>;
    static Type get_default(size_t size)
    {
        return Type(size, 50.0);
    }
    const static std::string name()
    {
        return "SeasonalitySigma";
    }
};

template <typename FP, class Status>
using ParametersBase = mio::ParameterSet<AdoptionRates<FP, Status>, TransitionRates<FP, Status>, FirstSeasonStartDay,
                                         StartDay, SeasonalityPeak, SeasonalityRho, SeasonalitySigma>;

} // namespace smm

} // namespace mio

#endif
