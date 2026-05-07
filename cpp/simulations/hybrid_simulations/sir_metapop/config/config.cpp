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
namespace sir
{

Config get_config(ConfigType type)
{
    double infected_commuting_reduc_fac;
    double transition_rate;
    Config config;
    switch (type) {
    case ConfigType::Config1:
        config.name              = "config1_1r_I0_1";
        config.num_regions       = 1;
        config.t0                = 0;
        config.tmax              = 90;
        config.dt                = 0.1;
        config.lambdas           = {0.00000007};
        config.gamma             = 1. / 7.;
        config.I0s               = {{0, 1}}; // Tuple of (region_id, initial infected)
        config.total_populations = {10000000};
        return config;
        break;
    case ConfigType::Config2:
        config.name              = "config1_1r_I0_2";
        config.num_regions       = 1;
        config.t0                = 0;
        config.tmax              = 90;
        config.dt                = 0.1;
        config.lambdas           = {0.00000007};
        config.gamma             = 1. / 7.;
        config.I0s               = {{0, 2}}; // Tuple of (region_id, initial infected)
        config.total_populations = {10000000};
        return config;
        break;
    case ConfigType::Config3:
        config.name              = "config1_1r_I0_10";
        config.num_regions       = 1;
        config.t0                = 0;
        config.tmax              = 90;
        config.dt                = 0.1;
        config.lambdas           = {0.00000007};
        config.gamma             = 1. / 7.;
        config.I0s               = {{0, 10}}; // Tuple of (region_id, initial infected)
        config.total_populations = {10000000};
        return config;
        break;
    case ConfigType::Config4:
        config.name              = "config1_1r_I0_100";
        config.num_regions       = 1;
        config.t0                = 0;
        config.tmax              = 90;
        config.dt                = 0.1;
        config.lambdas           = {0.00000007};
        config.gamma             = 1. / 7.;
        config.I0s               = {{0, 100}}; // Tuple of (region_id, initial infected)
        config.total_populations = {10000000};
        return config;
        break;
    case ConfigType::Config5:
        config.name              = "config1_1r_I0_1000";
        config.num_regions       = 1;
        config.t0                = 0;
        config.tmax              = 90;
        config.dt                = 0.1;
        config.lambdas           = {0.00000007};
        config.gamma             = 1. / 7.;
        config.I0s               = {{0, 1000}}; // Tuple of (region_id, initial infected)
        config.total_populations = {10000000};
        return config;
        break;
    case ConfigType::Config6:
        config.name              = "config1_1r_I0_5";
        config.num_regions       = 1;
        config.t0                = 0;
        config.tmax              = 90;
        config.dt                = 0.1;
        config.lambdas           = {0.00000007};
        config.gamma             = 1. / 7.;
        config.I0s               = {{0, 5}}; // Tuple of (region_id, initial infected)
        config.total_populations = {10000000};
        return config;
        break;
    case ConfigType::Config2regionsNoExchange:
        config.name                  = "config_2r_10_100_no_exchange";
        config.num_regions           = 2;
        config.t0                    = 0;
        config.tmax                  = 300;
        config.dt                    = 0.1;
        config.lambdas               = {0.00000007, 0.00000007};
        config.gamma                 = 1. / 7.;
        config.I0s                   = {{0, 100}, {1, 10}};
        config.total_populations     = {10000000, 10000000};
        infected_commuting_reduc_fac = 1.0;
        transition_rate              = 0;
        return config;
        break;
    case ConfigType::Config2regionsNoExchange2:
        config.name                  = "config_2r_5_100_no_exchange";
        config.num_regions           = 2;
        config.t0                    = 0;
        config.tmax                  = 300;
        config.dt                    = 0.1;
        config.lambdas               = {0.00000007, 0.00000007};
        config.gamma                 = 1. / 7.;
        config.I0s                   = {{0, 100}, {1, 5}};
        config.total_populations     = {10000000, 10000000};
        infected_commuting_reduc_fac = 1.0;
        transition_rate              = 0;
        return config;
        break;
    case ConfigType::Config2regionsk1:
        config.name                  = "config_2r_10_100_k1";
        config.num_regions           = 2;
        config.t0                    = 0;
        config.tmax                  = 90;
        config.dt                    = 0.1;
        config.lambdas               = {0.00000007, 0.00000007};
        config.gamma                 = 1. / 7.;
        config.I0s                   = {{0, 100}, {1, 10}};
        config.total_populations     = {10000000, 10000000};
        infected_commuting_reduc_fac = 1.0;
        transition_rate              = 0.0000001;
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
    case ConfigType::Config2regionsk2:
        config.name                  = "config_2r_10_100_k2";
        config.num_regions           = 2;
        config.t0                    = 0;
        config.tmax                  = 100;
        config.dt                    = 0.1;
        config.lambdas               = {0.00000007, 0.00000007};
        config.gamma                 = 1. / 7.;
        config.I0s                   = {{0, 100}, {1, 10}};
        config.total_populations     = {5000000, 5000000};
        infected_commuting_reduc_fac = 1.0;
        transition_rate              = 0.000001;
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
    case ConfigType::ConfigDiseaseImport1k1:
        config.name                  = "config_2r_100_0_k1";
        config.num_regions           = 2;
        config.t0                    = 0;
        config.tmax                  = 300;
        config.dt                    = 0.1;
        config.lambdas               = {0.00000007, 0.00000007};
        config.gamma                 = 1. / 7.;
        config.I0s                   = {{0, 100}, {1, 0}};
        config.total_populations     = {5000000, 5000000};
        infected_commuting_reduc_fac = 1.0;
        transition_rate              = 0.0000001;
        config.transition_rates.push_back({mio::osir::InfectionState::Infected, mio::regions::Region(0),
                                           mio::regions::Region(1), infected_commuting_reduc_fac * transition_rate});
        return config;
        break;
    case ConfigType::ConfigDiseaseImport1k2:
        config.name                  = "config_2r_100_0_k2";
        config.num_regions           = 2;
        config.t0                    = 0;
        config.tmax                  = 300;
        config.dt                    = 0.1;
        config.lambdas               = {0.00000007, 0.00000007};
        config.gamma                 = 1. / 7.;
        config.I0s                   = {{0, 100}, {1, 0}};
        config.total_populations     = {5000000, 5000000};
        infected_commuting_reduc_fac = 1.0;
        transition_rate              = 0.000001;
        config.transition_rates.push_back({mio::osir::InfectionState::Infected, mio::regions::Region(0),
                                           mio::regions::Region(1), infected_commuting_reduc_fac * transition_rate});
        return config;
        break;
    case ConfigType::ConfigDiseaseImport2k1:
        config.name                  = "config_2r_10_0_k1";
        config.num_regions           = 2;
        config.t0                    = 0;
        config.tmax                  = 300;
        config.dt                    = 0.1;
        config.lambdas               = {0.00000007, 0.00000007};
        config.gamma                 = 1. / 7.;
        config.I0s                   = {{0, 10}, {1, 0}};
        config.total_populations     = {5000000, 5000000};
        infected_commuting_reduc_fac = 1.0;
        transition_rate              = 0.0000001;
        config.transition_rates.push_back({mio::osir::InfectionState::Infected, mio::regions::Region(0),
                                           mio::regions::Region(1), infected_commuting_reduc_fac * transition_rate});
        return config;
        break;
    case ConfigType::ConfigDiseaseImport2k2:
        config.name                  = "config_2r_10_0_k2";
        config.num_regions           = 2;
        config.t0                    = 0;
        config.tmax                  = 300;
        config.dt                    = 0.1;
        config.lambdas               = {0.00000007, 0.00000007};
        config.gamma                 = 1. / 7.;
        config.I0s                   = {{0, 10}, {1, 0}};
        config.total_populations     = {5000000, 5000000};
        infected_commuting_reduc_fac = 1.0;
        transition_rate              = 0.000001;
        config.transition_rates.push_back({mio::osir::InfectionState::Infected, mio::regions::Region(0),
                                           mio::regions::Region(1), infected_commuting_reduc_fac * transition_rate});
        return config;
        break;
    case ConfigType::Config2regionsasymmExchange:
        config.name                  = "config_2r_100_10_k1_asym";
        config.num_regions           = 2;
        config.t0                    = 0;
        config.tmax                  = 100;
        config.dt                    = 0.1;
        config.lambdas               = {0.00000007, 0.00000007};
        config.gamma                 = 1. / 7.;
        config.I0s                   = {{0, 100}, {1, 10}};
        config.total_populations     = {5000000, 5000000};
        infected_commuting_reduc_fac = 1.0;
        transition_rate              = 0.0000001;
        config.transition_rates.push_back({mio::osir::InfectionState::Infected, mio::regions::Region(0),
                                           mio::regions::Region(1), infected_commuting_reduc_fac * transition_rate});
        return config;
        break;
    case ConfigType::Config4regionsk2:
        config.name                  = "config_4r_0_1_10_100_k2";
        config.num_regions           = 4;
        config.t0                    = 0;
        config.tmax                  = 300;
        config.dt                    = 0.1;
        config.lambdas               = {0.00000007, 0.00000007, 0.00000007, 0.00000007};
        config.gamma                 = 1. / 7.;
        config.I0s                   = {{0, 0}, {1, 1}, {2, 10}, {3, 100}};
        config.total_populations     = {5000000, 5000000, 5000000, 5000000};
        infected_commuting_reduc_fac = 1.0;
        transition_rate              = 0.000001;
        for (size_t region_from = 0; region_from < 4; ++region_from) {
            for (size_t region_to = 0; region_to < 4; ++region_to) {
                if (region_from != region_to) {
                    config.transition_rates.push_back({mio::osir::InfectionState::Susceptible,
                                                       mio::regions::Region(region_from),
                                                       mio::regions::Region(region_to), transition_rate});
                    config.transition_rates.push_back({mio::osir::InfectionState::Infected,
                                                       mio::regions::Region(region_from),
                                                       mio::regions::Region(region_to), transition_rate});
                    config.transition_rates.push_back({mio::osir::InfectionState::Recovered,
                                                       mio::regions::Region(region_from),
                                                       mio::regions::Region(region_to), transition_rate});
                }
            }
        }
        return config;
        break;
    case ConfigType::Config4regionsTransmNoExchange:
        config.name                  = "config_4r_10_transm_no_exchange";
        config.num_regions           = 4;
        config.t0                    = 0;
        config.tmax                  = 300;
        config.dt                    = 0.1;
        config.lambdas               = {0.00000004, 0.00000005, 0.00000008, 0.0000001};
        config.gamma                 = 1. / 7.;
        config.I0s                   = {{0, 10}, {1, 10}, {2, 10}, {3, 10}};
        config.total_populations     = {5000000, 5000000, 5000000, 5000000};
        infected_commuting_reduc_fac = 1.0;
        transition_rate              = 0.0;

        return config;
        break;
    case ConfigType::Config4regionsTransmSymExchange:
        config.name                  = "config_4r_10_transm_k2";
        config.num_regions           = 4;
        config.t0                    = 0;
        config.tmax                  = 300;
        config.dt                    = 0.1;
        config.lambdas               = {0.00000004, 0.00000005, 0.00000008, 0.0000001};
        config.gamma                 = 1. / 7.;
        config.I0s                   = {{0, 10}, {1, 10}, {2, 10}, {3, 10}};
        config.total_populations     = {5000000, 5000000, 5000000, 5000000};
        infected_commuting_reduc_fac = 1.0;
        transition_rate              = 0.000001;
        for (size_t region_from = 0; region_from < 4; ++region_from) {
            for (size_t region_to = 0; region_to < 4; ++region_to) {
                if (region_from != region_to) {
                    config.transition_rates.push_back({mio::osir::InfectionState::Susceptible,
                                                       mio::regions::Region(region_from),
                                                       mio::regions::Region(region_to), transition_rate});
                    config.transition_rates.push_back({mio::osir::InfectionState::Infected,
                                                       mio::regions::Region(region_from),
                                                       mio::regions::Region(region_to), transition_rate});
                    config.transition_rates.push_back({mio::osir::InfectionState::Recovered,
                                                       mio::regions::Region(region_from),
                                                       mio::regions::Region(region_to), transition_rate});
                }
            }
        }
        return config;
        break;
    case ConfigType::ConfigConference1:
        config.name        = "R0_4.9_I0_10";
        config.num_regions = 1;
        config.t0          = 0;
        config.tmax        = 100;
        config.dt          = 0.1;
        config.lambdas     = {0.00000007}; // 2 - 0.0000000286, 1.5 - 0.0000000215, 1 - 0.0000000143, 4.9 - 0.00000007
        config.gamma       = 1. / 7.;
        config.nu          = 1. / 20.;
        config.I0s         = {{0, 10}}; // Tuple of (region_id, initial infected)
        config.total_populations = {10000000};
        return config;
        break;
    }
    return config;
}
} // namespace sir

namespace sirs
{

Config get_config(ConfigType type)
{
    Config config;
    switch (type) {
    case ConfigType::Config1:
        config.name              = "config1_1r_I0_1";
        config.num_regions       = 1;
        config.t0                = 0;
        config.tmax              = 200;
        config.dt                = 0.1;
        config.lambdas           = {0.00007};
        config.gamma             = 1. / 7.;
        config.nu                = 1. / 20.;
        config.I0s               = {{0, 1}}; // Tuple of (region_id, initial infected)
        config.total_populations = {10000};
        return config;
        break;
    case ConfigType::Config3:
        config.name              = "config1_1r_I0_10";
        config.num_regions       = 1;
        config.t0                = 0;
        config.tmax              = 200;
        config.dt                = 0.1;
        config.lambdas           = {0.00007};
        config.gamma             = 1. / 7.;
        config.nu                = 1. / 20.;
        config.I0s               = {{0, 10}}; // Tuple of (region_id, initial infected)
        config.total_populations = {10000};
        return config;
        break;
    case ConfigType::Config4:
        config.name              = "config1_1r_I0_100";
        config.num_regions       = 1;
        config.t0                = 0;
        config.tmax              = 200;
        config.dt                = 0.1;
        config.lambdas           = {0.00007};
        config.gamma             = 1. / 7.;
        config.nu                = 1. / 20.;
        config.I0s               = {{0, 100}}; // Tuple of (region_id, initial infected)
        config.total_populations = {10000};
        return config;
        break;
    case ConfigType::Config4regionsTransmNoExchange:
        config.name              = "config_4r_10_transm_no_exchange";
        config.num_regions       = 4;
        config.t0                = 0;
        config.tmax              = 300;
        config.dt                = 0.1;
        config.lambdas           = {0.000021, 0.000028, 0.000042, 0.000057};
        config.gamma             = 1. / 7.;
        config.nu                = 1. / 20.;
        config.I0s               = {{0, 10}, {1, 10}, {2, 10}, {3, 10}};
        config.total_populations = {10000, 10000, 10000, 100000};

        return config;
        break;
    }
    return config;
}

} // namespace sirs

} // namespace Config
