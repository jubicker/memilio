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
#include "config.h"
#include <string>
#include <vector>

namespace Config
{

Config get_config(ConfigType type)
{
    Config config;
    switch (type) {
    case ConfigType::Config1:
        config.name              = "config1_1r_I0_1";
        config.num_regions       = 1;
        config.t0                = 0;
        config.tmax              = 90;
        config.dt                = 0.01;
        config.lambdas           = {0.00000007};
        config.gamma             = 1. / 7.;
        config.I0s               = {{0, 1}}; // Tuple of (region_id, initial infected)
        config.total_populations = {10000000};
        return config;
        break;
    case ConfigType::Config2:
        config.name              = "config1_1r_I0_10";
        config.num_regions       = 1;
        config.t0                = 0;
        config.tmax              = 90;
        config.dt                = 0.01;
        config.lambdas           = {0.00000007};
        config.gamma             = 1. / 7.;
        config.I0s               = {{0, 10}}; // Tuple of (region_id, initial infected)
        config.total_populations = {10000000};
        return config;
        break;
    case ConfigType::Config3:
        config.name              = "config1_1r_I0_100";
        config.num_regions       = 1;
        config.t0                = 0;
        config.tmax              = 90;
        config.dt                = 0.01;
        config.lambdas           = {0.00000007};
        config.gamma             = 1. / 7.;
        config.I0s               = {{0, 100}}; // Tuple of (region_id, initial infected)
        config.total_populations = {10000000};
        return config;
        break;
    case ConfigType::Config4:
        config.name              = "config2_1r_I0_1";
        config.num_regions       = 1;
        config.t0                = 0;
        config.tmax              = 90;
        config.dt                = 0.01;
        config.lambdas           = {0.00000014};
        config.gamma             = 1. / 7.;
        config.I0s               = {{0, 1}}; // Tuple of (region_id, initial infected)
        config.total_populations = {10000000};
        return config;
        break;
    case ConfigType::Config5:
        config.name              = "config2_1r_I0_10";
        config.num_regions       = 1;
        config.t0                = 0;
        config.tmax              = 90;
        config.dt                = 0.01;
        config.lambdas           = {0.00000014};
        config.gamma             = 1. / 7.;
        config.I0s               = {{0, 10}}; // Tuple of (region_id, initial infected)
        config.total_populations = {10000000};
        return config;
        break;
    case ConfigType::Config6:
        config.name              = "config2_1r_I0_100";
        config.num_regions       = 1;
        config.t0                = 0;
        config.tmax              = 90;
        config.dt                = 0.01;
        config.lambdas           = {0.00000014};
        config.gamma             = 1. / 7.;
        config.I0s               = {{0, 100}}; // Tuple of (region_id, initial infected)
        config.total_populations = {10000000};
        return config;
        break;
    case ConfigType::Config7:
        config.name              = "config1_1r_I0_2";
        config.num_regions       = 1;
        config.t0                = 0;
        config.tmax              = 90;
        config.dt                = 0.01;
        config.lambdas           = {0.00000007};
        config.gamma             = 1. / 7.;
        config.I0s               = {{0, 2}}; // Tuple of (region_id, initial infected)
        config.total_populations = {10000000};
        return config;
        break;
    case ConfigType::Config8:
        config.name              = "config2_1r_I0_2";
        config.num_regions       = 1;
        config.t0                = 0;
        config.tmax              = 90;
        config.dt                = 0.01;
        config.lambdas           = {0.00000014};
        config.gamma             = 1. / 7.;
        config.I0s               = {{0, 2}}; // Tuple of (region_id, initial infected)
        config.total_populations = {10000000};
        return config;
        break;
    case ConfigType::Config2r1:
        config.name                         = "config_2r1";
        config.num_regions                  = 2;
        config.t0                           = 0;
        config.tmax                         = 90;
        config.dt                           = 0.01;
        config.lambdas                      = {0.0000002, 0.0000002};
        config.gamma                        = 1. / 7.;
        config.I0s                          = {{0, 1}, {1, 0}};
        config.total_populations            = {5000000, 5000000};
        double infected_commuting_reduc_fac = 1.0;
        double transition_rate              = 0.00001;
        config.transition_rates.push_back({mio::osir::InfectionState::Susceptible, mio::regions::Region(0),
                                           mio::regions::Region(1), transition_rate});
        config.transition_rates.push_back({mio::osir::InfectionState::Susceptible, mio::regions::Region(1),
                                           mio::regions::Region(0), transition_rate});
        config.transition_rates.push_back(
            {mio::osir::InfectionState::Recovered, mio::regions::Region(0), mio::regions::Region(1), transition_rate});
        config.transition_rates.push_back(
            {mio::osir::InfectionState::Recovered, mio::regions::Region(1), mio::regions::Region(0), transition_rate});
        config.transition_rates.push_back({mio::osir::InfectionState::Infected, mio::regions::Region(0),
                                           mio::regions::Region(1), infected_commuting_reduc_fac * transition_rate});
        config.transition_rates.push_back({mio::osir::InfectionState::Infected, mio::regions::Region(1),
                                           mio::regions::Region(0), infected_commuting_reduc_fac * transition_rate});
        return config;
        break;
    }
    return config;
}

} // namespace Config
