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
#include <cstddef>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace Config
{

Config get_config(ConfigType type)
{
    double transition_rate;
    std::map<std::pair<size_t, size_t>, double> rates;
    Config config;
    switch (type) {
    case ConfigType::ConfigSIRVaryI0NoExchange:
        config.name                   = "config_SIR_I0_0_1_10_100_no_exchange";
        config.num_regions            = 4;
        config.t0                     = 0;
        config.tmax                   = 200;
        config.dt                     = 0.1;
        config.lambdas                = {0.00000286, 0.00000286, 0.00000286, 0.00000286}; //R0~2
        config.influencing_regions    = {{{0, 1}}, {{1, 1}}, {{2, 1}}, {{3, 1}}};
        config.gamma                  = 1. / 7.;
        config.nu                     = 0; // unused for SIR
        config.I0s                    = {{0, 0}, {1, 1}, {2, 10}, {3, 100}}; // Tuple of (region_id, initial infected)
        config.R0s                    = {}; // Tuple of (region_id, initial recovered)
        config.total_populations      = {100000, 100000, 100000, 100000};
        config.first_season_start_day = 0; // unused for SIR
        config.season_peaks           = {}; // unused for SIR
        config.seasonality_rhos       = {}; // unused for SIR
        config.seasonality_sigmas     = {}; // unused for SIR
        return config;
        break;
    case ConfigType::ConfigSIRVaryI0Exchange:
        config.name                   = "config_SIR_I0_0_1_10_100_exchange";
        config.num_regions            = 4;
        config.t0                     = 0;
        config.tmax                   = 200;
        config.dt                     = 0.1;
        config.lambdas                = {0.00000286, 0.00000286, 0.00000286, 0.00000286}; //R0~2
        config.influencing_regions    = {{{0, 1}}, {{1, 1}}, {{2, 1}}, {{3, 1}}};
        config.gamma                  = 1. / 7.;
        config.nu                     = 0; // unused for SIR
        config.I0s                    = {{0, 0}, {1, 1}, {2, 10}, {3, 100}}; // Tuple of (region_id, initial infected)
        config.R0s                    = {}; // Tuple of (region_id, initial recovered)
        config.total_populations      = {100000, 100000, 100000, 100000};
        config.first_season_start_day = 0; // unused for SIR
        config.season_peaks           = {}; // unused for SIR
        config.seasonality_rhos       = {}; // unused for SIR
        config.seasonality_sigmas     = {}; // unused for SIR
        transition_rate               = 0.00001;
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
    case ConfigType::ConfigSIRVaryR0NoExchange:
        config.name                   = "config_SIR_R0_1_1.5_2_4_no_exchange";
        config.num_regions            = 4;
        config.t0                     = 0;
        config.tmax                   = 200;
        config.dt                     = 0.1;
        config.lambdas                = {0.000001432, 0.00000215, 0.00000286, 0.00000572}; //R0~1, 1.5, 2, 4
        config.influencing_regions    = {{{0, 1}}, {{1, 1}}, {{2, 1}}, {{3, 1}}};
        config.gamma                  = 1. / 7.;
        config.nu                     = 0; // unused for SIR
        config.I0s                    = {{0, 10}, {1, 10}, {2, 10}, {3, 10}}; // Tuple of (region_id, initial infected)
        config.R0s                    = {}; // Tuple of (region_id, initial recovered)
        config.total_populations      = {100000, 100000, 100000, 100000};
        config.first_season_start_day = 0; // unused for SIR
        config.season_peaks           = {}; // unused for SIR
        config.seasonality_rhos       = {}; // unused for SIR
        config.seasonality_sigmas     = {}; // unused for SIR
        return config;
        break;
    case ConfigType::ConfigSIRVaryR0Exchange:
        config.name                   = "config_SIR_R0_1_1.5_2_4_exchange";
        config.num_regions            = 4;
        config.t0                     = 0;
        config.tmax                   = 200;
        config.dt                     = 0.1;
        config.lambdas                = {0.000001432, 0.00000215, 0.00000286, 0.00000572}; //R0~1, 1.5, 2, 4
        config.influencing_regions    = {{{0, 1}}, {{1, 1}}, {{2, 1}}, {{3, 1}}};
        config.gamma                  = 1. / 7.;
        config.nu                     = 0; // unused for SIR
        config.I0s                    = {{0, 10}, {1, 10}, {2, 10}, {3, 10}}; // Tuple of (region_id, initial infected)
        config.R0s                    = {}; // Tuple of (region_id, initial recovered)
        config.total_populations      = {100000, 100000, 100000, 100000};
        config.first_season_start_day = 0; // unused for SIR
        config.season_peaks           = {}; // unused for SIR
        config.seasonality_rhos       = {}; // unused for SIR
        config.seasonality_sigmas     = {}; // unused for SIR
        transition_rate               = 0.00001;
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
    case ConfigType::ConfigSIRSVaryI0NoExchange:
        config.name                   = "config_SIRS_I0_0_1_10_100_no_exchange";
        config.num_regions            = 4;
        config.t0                     = 0;
        config.tmax                   = 200;
        config.dt                     = 0.1;
        config.lambdas                = {0.00000286, 0.00000286, 0.00000286, 0.00000286}; //R0~2
        config.influencing_regions    = {{{0, 1}}, {{1, 1}}, {{2, 1}}, {{3, 1}}};
        config.gamma                  = 1. / 7.;
        config.nu                     = 1 / 20.;
        config.I0s                    = {{0, 0}, {1, 1}, {2, 10}, {3, 100}}; // Tuple of (region_id, initial infected)
        config.R0s                    = {}; // Tuple of (region_id, initial recovered)
        config.total_populations      = {100000, 100000, 100000, 100000};
        config.first_season_start_day = 0; // unused for SIRS
        config.season_peaks           = {}; // unused for SIRS
        config.seasonality_rhos       = {}; // unused for SIRS
        config.seasonality_sigmas     = {}; // unused for SIRS
        return config;
        break;
    case ConfigType::ConfigSIRSVaryI0Exchange:
        config.name                   = "config_SIRS_I0_0_1_10_100_exchange";
        config.num_regions            = 4;
        config.t0                     = 0;
        config.tmax                   = 200;
        config.dt                     = 0.1;
        config.lambdas                = {0.00000286, 0.00000286, 0.00000286, 0.00000286}; //R0~2
        config.influencing_regions    = {{{0, 1}}, {{1, 1}}, {{2, 1}}, {{3, 1}}};
        config.gamma                  = 1. / 7.;
        config.nu                     = 1 / 20.;
        config.I0s                    = {{0, 0}, {1, 1}, {2, 10}, {3, 100}}; // Tuple of (region_id, initial infected)
        config.R0s                    = {}; // Tuple of (region_id, initial recovered)
        config.total_populations      = {100000, 100000, 100000, 100000};
        config.first_season_start_day = 0; // unused for SIRS
        config.season_peaks           = {}; // unused for SIRS
        config.seasonality_rhos       = {}; // unused for SIRS
        config.seasonality_sigmas     = {}; // unused for SIRS
        transition_rate               = 0.00001;
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
    case ConfigType::ConfigSIRSVaryR0NoExchange:
        config.name                   = "config_SIRS_R0_1_1.5_2_4_no_exchange";
        config.num_regions            = 4;
        config.t0                     = 0;
        config.tmax                   = 200;
        config.dt                     = 0.1;
        config.lambdas                = {0.000001432, 0.00000215, 0.00000286, 0.00000572}; ////R0~1, 1.5, 2, 4
        config.influencing_regions    = {{{0, 1}}, {{1, 1}}, {{2, 1}}, {{3, 1}}};
        config.gamma                  = 1. / 7.;
        config.nu                     = 1 / 20.;
        config.I0s                    = {{0, 10}, {1, 10}, {2, 10}, {3, 10}}; // Tuple of (region_id, initial infected)
        config.R0s                    = {}; // Tuple of (region_id, initial recovered)
        config.total_populations      = {100000, 100000, 100000, 100000};
        config.first_season_start_day = 0; // unused for SIRS
        config.season_peaks           = {}; // unused for SIRS
        config.seasonality_rhos       = {}; // unused for SIRS
        config.seasonality_sigmas     = {}; // unused for SIRS
        return config;
        break;
    case ConfigType::ConfigSIRSVaryR0Exchange:
        config.name                   = "config_SIRS_R0_1_1.5_2_4_exchange";
        config.num_regions            = 4;
        config.t0                     = 0;
        config.tmax                   = 200;
        config.dt                     = 0.1;
        config.lambdas                = {0.000001432, 0.00000215, 0.00000286, 0.00000572}; //R0~1, 1.5, 2, 4
        config.influencing_regions    = {{{0, 1}}, {{1, 1}}, {{2, 1}}, {{3, 1}}};
        config.gamma                  = 1. / 7.;
        config.nu                     = 1 / 20.;
        config.I0s                    = {{0, 10}, {1, 10}, {2, 10}, {3, 10}}; // Tuple of (region_id, initial infected)
        config.R0s                    = {}; // Tuple of (region_id, initial recovered)
        config.total_populations      = {100000, 100000, 100000, 100000};
        config.first_season_start_day = 0; // unused for SIRS
        config.season_peaks           = {}; // unused for SIRS
        config.seasonality_rhos       = {}; // unused for SIRS
        config.seasonality_sigmas     = {}; // unused for SIRS
        transition_rate               = 0.00001;
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
    case ConfigType::ConfigSIRSVaryI0NoExchange_Seasonal:
        config.name                   = "config_SIRS_seasonal_I0_0_1_10_100_no_exchange";
        config.num_regions            = 4;
        config.t0                     = 0;
        config.tmax                   = 3 * 365; // 3 years to see multiple seasons
        config.dt                     = 0.1;
        config.lambdas                = {0.00000572, 0.00000572, 0.00000572, 0.00000572}; //R0~4
        config.influencing_regions    = {{{0, 1}}, {{1, 1}}, {{2, 1}}, {{3, 1}}};
        config.gamma                  = 1. / 7.;
        config.nu                     = 1 / 20.;
        config.I0s                    = {{0, 0}, {1, 1}, {2, 10}, {3, 100}}; // Tuple of (region_id, initial infected)
        config.R0s                    = {}; // Tuple of (region_id, initial recovered)
        config.total_populations      = {100000, 100000, 100000, 100000};
        config.first_season_start_day = -180; // First season starts at 1st July of the previous year
        config.season_peaks           = {24, 24, 24}; // Peak transmission rate at day January 24th of each season
        config.seasonality_rhos       = {0.4, 0.4, 0.4}; // Variation factor for each season
        config.seasonality_sigmas     = {50, 50, 50}; // Standard deviation for the Gaussian function for each season
        return config;
        break;
    case ConfigType::ConfigSIRSVaryI0Exchange_Seasonal:
        config.name                   = "config_SIRS_seasonal_I0_0_1_10_100_exchange";
        config.num_regions            = 4;
        config.t0                     = 0;
        config.tmax                   = 3 * 365; // 3 years to see multiple seasons
        config.dt                     = 0.1;
        config.lambdas                = {0.00000572, 0.00000572, 0.00000572, 0.00000572}; //R0~4
        config.influencing_regions    = {{{0, 1}}, {{1, 1}}, {{2, 1}}, {{3, 1}}};
        config.gamma                  = 1. / 7.;
        config.nu                     = 1 / 20.;
        config.I0s                    = {{0, 0}, {1, 1}, {2, 10}, {3, 100}}; // Tuple of (region_id, initial infected)
        config.R0s                    = {}; // Tuple of (region_id, initial recovered)
        config.total_populations      = {100000, 100000, 100000, 100000};
        config.first_season_start_day = -180; // First season starts at 1st July of the previous year
        config.season_peaks           = {24, 24, 24}; // Peak transmission rate at day January 24th of each season
        config.seasonality_rhos       = {0.4, 0.4, 0.4}; // Variation factor for each season
        config.seasonality_sigmas     = {50, 50, 50}; // Standard deviation for the Gaussian function for each season
        transition_rate               = 0.00001;
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
    case ConfigType::ConfigSIRSVaryR0NoExchange_Seasonal:
        config.name                   = "config_SIRS_seasonal_R0_1_1.5_2_4_no_exchange";
        config.num_regions            = 4;
        config.t0                     = 0;
        config.tmax                   = 3 * 365;
        config.dt                     = 0.1;
        config.lambdas                = {0.000001432, 0.00000215, 0.00000286, 0.00000572}; ////R0~1, 1.5, 2, 4
        config.influencing_regions    = {{{0, 1}}, {{1, 1}}, {{2, 1}}, {{3, 1}}};
        config.gamma                  = 1. / 7.;
        config.nu                     = 1 / 20.;
        config.I0s                    = {{0, 10}, {1, 10}, {2, 10}, {3, 10}}; // Tuple of (region_id, initial infected)
        config.R0s                    = {}; // Tuple of (region_id, initial recovered)
        config.total_populations      = {100000, 100000, 100000, 100000};
        config.first_season_start_day = -180; // First season starts at 1st July of the previous year
        config.season_peaks           = {24, 24, 24}; // Peak transmission rate at day January 24th of each season
        config.seasonality_rhos       = {0.4, 0.4, 0.4}; // Variation factor for each season
        config.seasonality_sigmas     = {50, 50, 50}; // Standard deviation for the Gaussian function for each season
        return config;
        break;
    case ConfigType::ConfigSIRSVaryR0Exchange_Seasonal:
        config.name                   = "config_SIRS_seasonal_R0_1_1.5_2_4_exchange";
        config.num_regions            = 4;
        config.t0                     = 0;
        config.tmax                   = 3 * 365;
        config.dt                     = 0.1;
        config.lambdas                = {0.000001432, 0.00000215, 0.00000286, 0.00000572}; ////R0~1, 1.5, 2, 4
        config.influencing_regions    = {{{0, 1}}, {{1, 1}}, {{2, 1}}, {{3, 1}}};
        config.gamma                  = 1. / 7.;
        config.nu                     = 1 / 20.;
        config.I0s                    = {{0, 10}, {1, 10}, {2, 10}, {3, 10}}; // Tuple of (region_id, initial infected)
        config.R0s                    = {}; // Tuple of (region_id, initial recovered)
        config.total_populations      = {100000, 100000, 100000, 100000};
        config.first_season_start_day = -180; // First season starts at 1st July of the previous year
        config.season_peaks           = {24, 24, 24}; // Peak transmission rate at day January 24th of each season
        config.seasonality_rhos       = {0.4, 0.4, 0.4}; // Variation factor for each season
        config.seasonality_sigmas     = {50, 50, 50}; // Standard deviation for the Gaussian function for each season
        transition_rate               = 0.00001;
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
    case ConfigType::ConfigTestExchange:
        config.name                   = "config_SIR_test_exchange";
        config.num_regions            = 2;
        config.t0                     = 0;
        config.tmax                   = 200;
        config.dt                     = 0.1;
        config.lambdas                = {0.00000572, 0.00000572}; //R0~4
        config.influencing_regions    = {{{0, 1}}, {{1, 1}}};
        config.gamma                  = 1. / 7.;
        config.nu                     = 0; // unused for SIR
        config.I0s                    = {{0, 10}, {1, 0}}; // Tuple of (region_id, initial infected)
        config.R0s                    = {}; // Tuple of (region_id, initial recovered)
        config.total_populations      = {100000, 100000};
        config.first_season_start_day = 0; // unused for SIR
        config.season_peaks           = {}; // unused for SIR
        config.seasonality_rhos       = {}; // unused for SIR
        config.seasonality_sigmas     = {}; // unused for SIR
        transition_rate               = 0.00001;

        config.transition_rates.push_back({mio::osir::InfectionState::Susceptible, mio::regions::Region(0),
                                           mio::regions::Region(1), transition_rate});
        config.transition_rates.push_back(
            {mio::osir::InfectionState::Infected, mio::regions::Region(0), mio::regions::Region(1), transition_rate});
        config.transition_rates.push_back(
            {mio::osir::InfectionState::Recovered, mio::regions::Region(0), mio::regions::Region(1), transition_rate});

        return config;
        break;
    case ConfigType::ConfigTest2ndOutbreak:
        config.name                   = "config_SIR_test_2nd_outbreak";
        config.num_regions            = 1;
        config.t0                     = 0;
        config.tmax                   = 200;
        config.dt                     = 0.1;
        config.lambdas                = {0.00000572}; //R0~4
        config.influencing_regions    = {{{0, 1}}};
        config.gamma                  = 1. / 7.;
        config.nu                     = 0; // unused for SIR
        config.I0s                    = {{0, 1}}; // Tuple of (region_id, initial infected)
        config.R0s                    = {}; // Tuple of (region_id, initial recovered)
        config.total_populations      = {100000};
        config.first_season_start_day = 0; // unused for SIR
        config.season_peaks           = {}; // unused for SIR
        config.seasonality_rhos       = {}; // unused for SIR
        config.seasonality_sigmas     = {}; // unused for SIR

        return config;
        break;
    case ConfigType::ConfigPerformanceStudySIR:
        config.name                   = "performance_study_SIR";
        config.num_regions            = 1;
        config.t0                     = 0;
        config.tmax                   = 200;
        config.dt                     = 0.1;
        config.lambdas                = {0};
        config.influencing_regions    = {{{0, 1}}};
        config.gamma                  = 1. / 7.;
        config.nu                     = 0; // unused for SIR
        config.I0s                    = {{0, 0}}; // Tuple of (region_id, initial infected)
        config.R0s                    = {{0, 0}}; // Tuple of (region_id, initial recovered)
        config.total_populations      = {100000};
        config.first_season_start_day = 0; // unused for SIR
        config.season_peaks           = {}; // unused for SIR
        config.seasonality_rhos       = {}; // unused for SIR
        config.seasonality_sigmas     = {}; // unused for SIR

        return config;
        break;
    case ConfigType::ConfigPerformanceStudySIRS:
        config.name                = "performance_study_SIRS";
        config.num_regions         = 1;
        config.t0                  = 0;
        config.tmax                = 200;
        config.dt                  = 0.1;
        config.lambdas             = {0};
        config.influencing_regions = {{{0, 1}}};
        config.gamma               = 1. / 7.;
        config.nu                  = 1 / 20.;
        ;
        config.I0s                    = {{0, 0}}; // Tuple of (region_id, initial infected)
        config.R0s                    = {{0, 0}}; // Tuple of (region_id, initial recovered)
        config.total_populations      = {100000};
        config.first_season_start_day = 0; // unused for SIRS
        config.season_peaks           = {}; // unused for SIRS
        config.seasonality_rhos       = {}; // unused for SIRS
        config.seasonality_sigmas     = {}; // unused for SIRS

        return config;
        break;
    case ConfigType::ConfigPerformanceStudySIRS_seasonal:
        config.name                = "performance_study_SIRS_seasonal";
        config.num_regions         = 1;
        config.t0                  = 0;
        config.tmax                = 200;
        config.dt                  = 0.1;
        config.lambdas             = {0};
        config.influencing_regions = {{{0, 1}}};
        config.gamma               = 1. / 7.;
        config.nu                  = 1 / 20.;
        ;
        config.I0s                    = {{0, 0}}; // Tuple of (region_id, initial infected)
        config.R0s                    = {{0, 0}}; // Tuple of (region_id, initial recovered)
        config.total_populations      = {100000};
        config.first_season_start_day = -180; // First season starts at 1st July of the previous year
        config.season_peaks           = {24, 24, 24}; // Peak transmission rate at day January 24th of each season
        config.seasonality_rhos       = {0.4, 0.4, 0.4}; // Variation factor for each season
        config.seasonality_sigmas     = {50, 50, 50}; // Standard deviation for the Gaussian function for each season

        return config;
        break;
    case ConfigType::ConfigInfluenzaGermany:
        config.name                = "config_influenza_germany";
        config.num_regions         = 1;
        config.t0                  = 0;
        config.tmax                = 4 * 365 - 137; // From 1st August 2016 to 27th July 2020
        config.dt                  = 0.1;
        config.lambdas             = {0.0000025145241}; // manual: 0.0000018, top: 0.0000025145241, mean: 0.000002390025
        config.influencing_regions = {{{0, 1}}};
        config.gamma               = 1. / 7.;
        config.nu                  = 1 / 149.; // manual: 1/14, top: 1/149, mean: 1/128
        config.I0s                 = {{0, 513}}; // manual: 1000, top: 513, mean: 480
        config.R0s                 = {{0, 23718}}; // manual: 20000, top: 23718, mean: 16920
        config.total_populations   = {100000};
        config.first_season_start_day = -153; // First season starts at 1st August of the previous year
        config.season_peaks = {21, 37, 24, 24, 24}; // manual: {-10, -25, -20}, top: {21, 37, 24}, mean: {22, 38, 26}
        config.seasonality_rhos = {
            0.79, 0.78, 0.77, 0.76,
            0.775}; // manual: {0.5, 0.65, 0.1}, top: {0.79, 0.78, 0.77}, mean: {0.78, 0.78, 0.77}
        config.seasonality_sigmas = {10, 15, 5, 7, 8.5}; // manual: {100, 105, 130}, top: {10, 15, 5}, mean: {9, 15, 6}
        return config;
        break;
    case ConfigType::ConfigInfluenzaGermany_full:
        config.name                   = "config_influenza_germany_full";
        config.num_regions            = 1;
        config.t0                     = 0;
        config.tmax                   = 5 * 365 - 4;
        config.dt                     = 0.1;
        config.lambdas                = {2.908628081137058e-7};
        config.influencing_regions    = {{{0, 1}}};
        config.gamma                  = 1. / 7.;
        config.nu                     = 1 / 149.;
        config.I0s                    = {{0, int(4407.50)}};
        config.R0s                    = {{0, int(198228.26)}};
        config.total_populations      = {int(835771.40)};
        config.first_season_start_day = -153;
        config.season_peaks           = {21, 37, 24, 24, 24};
        config.seasonality_rhos       = {0.81, 0.83, 0.8, 0.78, 0.805};
        config.seasonality_sigmas     = {20, 15, 10, 10, 12.5};
        return config;
        break;
    case ConfigType::ConfigInfluenzaAgeGroups:
        config.name        = "config_influenza_agegroups";
        config.num_regions = 5; // 0-4, 5-14, 15-34, 35-59, 60+
        config.t0          = 0;
        config.tmax        = 3 * 365 - 4;
        config.dt          = 0.1;
        config.lambdas     = {0.00001288089301181601, 1.19954495503683e-6, 1.115361116626335e-6, 8.44593364010281e-7,
                              9.793771309493856e-7};
        config.influencing_regions    = {{{0, 0.0310993759369821},
                                          {1, 0.3220221286967569},
                                          {2, 0.6114921106569616},
                                          {3, 0.617696650073716},
                                          {4, 0.2638556578331692}},
                                         {{0, 0.0154584335258646},
                                          {1, 0.8357863483840422},
                                          {2, 0.7343745784912392},
                                          {3, 0.6873489020412964},
                                          {4, 0.5982126107792187}},
                                         {{0, 0.0035934065966545},
                                          {1, 0.5790198909973192},
                                          {2, 0.3920762942506107},
                                          {3, 0.6179550746528993},
                                          {4, 0.4955227268650907}},
                                         {{0, 0.004435413822006},
                                          {1, 0.7726787052175906},
                                          {2, 0.6331439884923853},
                                          {3, 0.1787871207107431},
                                          {4, 0.5657142406035242}},
                                         {{0, 0.0044687206750555},
                                          {1, 0.1949804086187295},
                                          {2, 0.9156397710855608},
                                          {3, 0.4809696098819842},
                                          {4, 0.0777360459261227}}};
        config.gamma                  = 1. / 7.;
        config.nu                     = 1 / 149.;
        config.I0s                    = {{0, 2181}, {1, 471}, {2, 98}, {3, 98}, {4, 177}};
        config.R0s                    = {{0, 9268}, {1, 29065}, {2, 25507}, {3, 29961}, {4, 13832}};
        config.total_populations      = {100000, 100000, 100000, 100000, 100000};
        config.first_season_start_day = -153;
        config.season_peaks           = {21, 37, 24};
        config.seasonality_rhos       = {0.79, 0.78, 0.77};
        config.seasonality_sigmas     = {10, 15, 5};
        return config;
        break;
    case ConfigType::ConfigInfluenzaAgeGroups_full:
        config.name                   = "config_influenza_agegroups_full";
        config.num_regions            = 5; // 0-4, 5-14, 15-34, 35-59, 60+
        config.t0                     = 0;
        config.tmax                   = 3 * 365 - 3;
        config.dt                     = 0.1;
        config.lambdas                = {2.04628894e-8, 1.39171521e-8, 3.7836113e-9, 3.40942323e-9, 1.9126576e-9};
        config.influencing_regions    = {{{0, 0.0016821903018872},
                                          {1, 0.0835973041131203},
                                          {2, 0.0594331854564333},
                                          {3, 0.6255548648395776},
                                          {4, 0.8325660776336196}},
                                         {{0, 0.0835973041131203},
                                          {1, 0.1090796848279727},
                                          {2, 0.5430772129054124},
                                          {3, 0.5387381290770932},
                                          {4, 0.3142973941728622}},
                                         {{0, 0.0594331854564333},
                                          {1, 0.5430772129054124},
                                          {2, 0.7980979140492013},
                                          {3, 0.918715326407079},
                                          {4, 0.8406336169560021}},
                                         {{0, 0.6255548648395776},
                                          {1, 0.5387381290770932},
                                          {2, 0.918715326407079},
                                          {3, 0.2952636345929488},
                                          {4, 0.78904291699207}},
                                         {{0, 0.8325660776336196},
                                          {1, 0.3142973941728622},
                                          {2, 0.8406336169560021},
                                          {3, 0.78904291699207},
                                          {4, 0.5242714373201034}}};
        config.gamma                  = 1. / 7.;
        config.nu                     = 1 / 149.;
        config.I0s                    = {{0, 58284}, {1, 45615}, {2, 33035}, {3, 146564}, {4, 42435}};
        config.R0s                    = {{0, 987247}, {1, 2226587}, {2, 4365956}, {3, 7161450}, {4, 5851266}};
        config.total_populations      = {3700607, 7919288, 18770235, 28077512, 25109498};
        config.first_season_start_day = -53;
        config.season_peaks           = {21, 37, 24};
        config.seasonality_rhos       = {0.79, 0.78, 0.77};
        config.seasonality_sigmas     = {10, 15, 5};
        return config;
        break;
    case ConfigType::ConfigInfluenzaRegions:
        config.name                   = "config_influenza_regions_regional_no_64cores";
        config.num_regions            = 4; // Norden (West), Osten, Sueden, Mitte (West)
        config.t0                     = 0;
        config.tmax                   = 5 * 365 - 4;
        config.dt                     = 0.1;
        config.lambdas                = {4.17489e-6, 3.161967e-6, 2.1127735e-6, 1.6121749e-6};
        config.influencing_regions    = {{{0, 1}}, {{1, 1}}, {{2, 1}}, {{3, 1}}};
        config.gamma                  = 1. / 7.;
        config.nu                     = 1 / 149.;
        config.I0s                    = {{0, 26}, {1, 60}, {2, 649}, {3, 763}};
        config.R0s                    = {{0, 28602}, {1, 29626}, {2, 19516}, {3, 19166}};
        config.total_populations      = {100000, 100000, 100000, 100000};
        config.first_season_start_day = -153;
        config.season_peaks           = {21, 37, 24, 24, 24};
        config.seasonality_rhos       = {0.79, 0.78, 0.77, 0.76, 0.775};
        config.seasonality_sigmas     = {10, 15, 5, 7, 8.5};
        // Fill transition rates - work counts
        rates = {{{0, 1}, 0.006735419081411219},  {{0, 3}, 0.011836497664847793},  {{1, 0}, 0.0056630134202979325},
                 {{1, 2}, 0.0028037943556893123}, {{1, 3}, 0.001536233748117688},  {{2, 1}, 0.0018421849577539354},
                 {{2, 3}, 0.010303359574793468},  {{3, 0}, 0.0054372554503847766}, {{3, 1}, 0.0008393263431793039},
                 {{3, 2}, 0.008567721370540751}};
        for (auto const& [key, val] : rates) {
            config.transition_rates.push_back({mio::osir::InfectionState::Susceptible, mio::regions::Region(key.first),
                                               mio::regions::Region(key.second), val});
            config.transition_rates.push_back({mio::osir::InfectionState::Infected, mio::regions::Region(key.first),
                                               mio::regions::Region(key.second), val});
            config.transition_rates.push_back({mio::osir::InfectionState::Recovered, mio::regions::Region(key.first),
                                               mio::regions::Region(key.second), val});
        }
        return config;
        break;
    case ConfigType::ConfigInfluenzaRegions_full:
        config.name        = "config_influenza_regions_full";
        config.num_regions = 4; // Norden (West), Osten, Sueden, Mitte (West)
        config.t0          = 0;
        config.tmax        = 5 * 365 - 4;
        config.dt          = 0.1;
        config.lambdas = {0.00000179924907097591, 0.000001531772435903319, 9.7867373864923955e-7, 8.152123248521465e-7};
        config.influencing_regions    = {{{0, 1}}, {{1, 1}}, {{2, 1}}, {{3, 1}}};
        config.gamma                  = 1. / 7.;
        config.nu                     = 1 / 149.;
        config.I0s                    = {{0, int(320)}, {1, int(1000)}, {2, int(1400)}, {3, int(1000)}};
        config.R0s                    = {{0, int(50000)}, {1, int(40220)}, {2, int(41850)}, {3, int(37700)}};
        config.total_populations      = {int(135314.52), int(160939.05), int(244948.26), int(294569.57)};
        config.first_season_start_day = -153;
        config.season_peaks           = {21, 37, 24, 24, 24};
        config.seasonality_rhos       = {0.81, 0.83, 0.8, 0.78, 0.805};
        config.seasonality_sigmas     = {20, 15, 10, 10, 12.5};
        // Fill transition rates - work counts
        rates = {{{0, 1}, 0.006735419081411219},  {{0, 3}, 0.011836497664847793},  {{1, 0}, 0.0056630134202979325},
                 {{1, 2}, 0.0028037943556893123}, {{1, 3}, 0.001536233748117688},  {{2, 1}, 0.0018421849577539354},
                 {{2, 3}, 0.010303359574793468},  {{3, 0}, 0.0054372554503847766}, {{3, 1}, 0.0008393263431793039},
                 {{3, 2}, 0.008567721370540751}};
        for (auto const& [key, val] : rates) {
            config.transition_rates.push_back({mio::osir::InfectionState::Susceptible, mio::regions::Region(key.first),
                                               mio::regions::Region(key.second), val});
            config.transition_rates.push_back({mio::osir::InfectionState::Infected, mio::regions::Region(key.first),
                                               mio::regions::Region(key.second), val});
            config.transition_rates.push_back({mio::osir::InfectionState::Recovered, mio::regions::Region(key.first),
                                               mio::regions::Region(key.second), val});
        }
        return config;
        break;
    case ConfigType::ConfigTest:
        config.name                = "eight_regions";
        config.num_regions         = 8;
        config.t0                  = 0;
        config.tmax                = 200;
        config.dt                  = 0.1;
        config.lambdas             = {0.000003608, 0.000003608, 0.000003608, 0.000003608,
                                      0.000003608, 0.000003608, 0.000003608, 0.000003608}; //R0~2.5
        config.influencing_regions = {{{0, 1}}, {{1, 1}}, {{2, 1}}, {{3, 1}}, {{4, 1}}, {{5, 1}}, {{6, 1}}, {{7, 1}}};
        config.gamma               = 1. / 7.;
        config.nu                  = 1 / 20.;
        config.I0s                 = {{0, 10}, {1, 10}, {2, 10}, {3, 10}, {4, 10}, {5, 10}, {6, 10}, {7, 10}};
        config.total_populations   = {100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000};
        config.first_season_start_day = 0; //-182;
        config.season_peaks           = {}; //{0, 0, 0};
        config.seasonality_rhos       = {}; //{0.5, 0.5, 0.5};
        config.seasonality_sigmas     = {}; //{40, 40, 40};
        transition_rate               = 0.00001;
        for (size_t region_from = 0; region_from < config.num_regions; ++region_from) {
            for (size_t region_to = 0; region_to < config.num_regions; ++region_to) {
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
    }
    return config;
}

} // namespace Config
