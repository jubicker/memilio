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
    double transition_rate;
    Config config;
    switch (type) {
    case ConfigType::ConfigSIRVaryI0NoExchange:
        config.name                   = "config_SIR_I0_0_1_10_100_no_exchange";
        config.num_regions            = 4;
        config.t0                     = 0;
        config.tmax                   = 200;
        config.dt                     = 0.1;
        config.lambdas                = {0.00000572, 0.00000572, 0.00000572, 0.00000572}; //R0~4
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
        config.lambdas                = {0.00000572, 0.00000572, 0.00000572, 0.00000572}; //R0~4
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
        config.lambdas                = {0.00000572, 0.00000572, 0.00000572, 0.00000572}; //R0~4
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
        config.lambdas                = {0.00000572, 0.00000572, 0.00000572, 0.00000572}; //R0~4
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
        config.name                   = "config_influenza_germany";
        config.num_regions            = 1;
        config.t0                     = 0;
        config.tmax                   = 3 * 365 - 3; // From 1st August 2016 to 29th July 2019
        config.dt                     = 0.1;
        config.lambdas                = {0.000003}; // TODO
        config.influencing_regions    = {{{0, 1}}};
        config.gamma                  = 1. / 7.;
        config.nu                     = 1 / 14.;
        config.I0s                    = {{0, 2500}}; // TODO
        config.R0s                    = {{0, 20000}}; // TODO
        config.total_populations      = {100000};
        config.first_season_start_day = -153; // First season starts at 1st August of the previous year
        config.season_peaks           = {22, 20, 20};
        config.seasonality_rhos       = {0.6, 0.6, 0.6}; // Variation factor for each season
        config.seasonality_sigmas     = {50, 50, 50}; // Standard deviation for the Gaussian function for each season
        return config;
        break;
    case ConfigType::ConfigInfluenzaGermany_full:
        config.name                   = "config_influenza_germany_full";
        config.num_regions            = 1;
        config.t0                     = 0;
        config.tmax                   = 3 * 365 + 4; // From 3rd October 2016 to 5th October 2019
        config.dt                     = 0.1;
        config.lambdas                = {0.000003}; // TODO
        config.influencing_regions    = {{{0, 1}}};
        config.gamma                  = 1. / 7.;
        config.nu                     = 1 / 14.;
        config.I0s                    = {{0, 1000}}; // TODO
        config.R0s                    = {{0, 20000}}; // TODO
        config.total_populations      = {83577140};
        config.first_season_start_day = -90; // First season starts at 3rd October of the previous year
        config.season_peaks           = {20, 20, 20};
        config.seasonality_rhos       = {0.3, 0.2, 0.2}; // Variation factor for each season
        config.seasonality_sigmas     = {40, 40, 40}; // Standard deviation for the Gaussian function for each season
        return config;
        break;
    case ConfigType::ConfigInfluenzaAgeGroups:
        config.name                   = "config_influenza_agegroups";
        config.num_regions            = 5; // 0-4, 5-14, 15-34, 35-59, 60+
        config.t0                     = 0;
        config.tmax                   = 3 * 365 + 4; // From 3rd October 2016 to 5th October 2019
        config.dt                     = 0.1;
        config.lambdas                = {0.000003, 0.000003, 0.000003, 0.000003, 0.000003}; // TODO
        config.influencing_regions    = {{{0, 1}, {1, 1}, {2, 1}, {3, 1}, {4, 1}},
                                         {{0, 1}, {1, 1}, {2, 1}, {3, 1}, {4, 1}},
                                         {{0, 1}, {1, 1}, {2, 1}, {3, 1}, {4, 1}},
                                         {{0, 1}, {1, 1}, {2, 1}, {3, 1}, {4, 1}},
                                         {{0, 1}, {1, 1}, {2, 1}, {3, 1}, {4, 1}}};
        config.gamma                  = 1. / 7.;
        config.nu                     = 1 / 14.;
        config.I0s                    = {{0, 1000}, {1, 1000}, {2, 1000}, {3, 1000}, {4, 1000}}; // TODO
        config.R0s                    = {{0, 20000}, {1, 20000}, {2, 20000}, {3, 20000}, {4, 20000}}; // TODO
        config.total_populations      = {100000, 100000, 100000, 100000, 100000};
        config.first_season_start_day = -90; // First season starts at 3rd October of the previous year
        config.season_peaks           = {20, 20, 20};
        config.seasonality_rhos       = {0.3, 0.2, 0.2}; // Variation factor for each season
        config.seasonality_sigmas     = {40, 40, 40}; // Standard deviation for the Gaussian function for each season
        return config;
        break;
    case ConfigType::ConfigInfluenzaAgeGroups_full:
        config.name                   = "config_influenza_agegroups_full";
        config.num_regions            = 5; // 0-4, 5-14, 15-34, 35-59, 60+
        config.t0                     = 0;
        config.tmax                   = 3 * 365 + 4; // From 3rd October 2016 to 5th October 2019
        config.dt                     = 0.1;
        config.lambdas                = {0.000003, 0.000003, 0.000003, 0.000003, 0.000003}; // TODO
        config.influencing_regions    = {{{0, 1}, {1, 1}, {2, 1}, {3, 1}, {4, 1}},
                                         {{0, 1}, {1, 1}, {2, 1}, {3, 1}, {4, 1}},
                                         {{0, 1}, {1, 1}, {2, 1}, {3, 1}, {4, 1}},
                                         {{0, 1}, {1, 1}, {2, 1}, {3, 1}, {4, 1}},
                                         {{0, 1}, {1, 1}, {2, 1}, {3, 1}, {4, 1}}};
        config.gamma                  = 1. / 7.;
        config.nu                     = 1 / 14.;
        config.I0s                    = {{0, 1000}, {1, 1000}, {2, 1000}, {3, 1000}, {4, 1000}}; // TODO
        config.R0s                    = {{0, 20000}, {1, 20000}, {2, 20000}, {3, 20000}, {4, 20000}}; // TODO
        config.total_populations      = {3700607, 7919288, 18770235, 28077512, 25109498};
        config.first_season_start_day = -90; // First season starts at 3rd October of the previous year
        config.season_peaks           = {20, 20, 20};
        config.seasonality_rhos       = {0.3, 0.2, 0.2}; // Variation factor for each season
        config.seasonality_sigmas     = {40, 40, 40}; // Standard deviation for the Gaussian function for each season
        return config;
        break;
    }
    return config;
}

} // namespace Config
