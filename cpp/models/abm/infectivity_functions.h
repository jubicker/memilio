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

#ifndef MIO_ABM_INFECTIVITY_FUNCTIONS_H
#define MIO_ABM_INFECTIVITY_FUNCTIONS_H

#include "abm/time.h"
#include "abm/infection.h"
#include "memilio/config.h"

namespace mio
{
namespace abm
{

/**
 * @brief Computes the infectivity at time t as sigmoid function of the viral load.
 * @param[in] t TimePoint of the query.
 * @param[in] infection Infection of which the infectivity should be calculated.
 */
ScalarType sigmoidal_infectivity(TimePoint t, const Infection& infection);

} // namespace abm
} // namespace mio

#endif // MIO_ABM_INFECTIVITY_FUNCTIONS_H
