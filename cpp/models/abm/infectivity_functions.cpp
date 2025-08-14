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

#include "abm/infectivity_functions.h"
#include "abm/infection.h"

namespace mio
{
namespace abm
{

ScalarType sigmoidal_infectivity(TimePoint t, const Infection& infection)
{
    return infection.get_virus_shed_factor() /
           (1 + exp(-(infection.get_alpha() + infection.get_beta() * infection.get_viral_load(t))));
}

ScalarType constant_infectivity(TimePoint /*t*/, const Infection& infection)
{
    return infection.get_virus_shed_factor();
}

} // namespace abm
} // namespace mio
