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
        config.name        = "performance_study_SIRS";
        config.num_regions = 1;
        config.t0          = 0;
        config.tmax        = 200;
        config.dt          = 0.1;
        config.lambdas     = {0};
        config.gamma       = 1. / 7.;
        config.nu          = 1 / 20.;
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
        config.name        = "performance_study_SIRS_seasonal";
        config.num_regions = 1;
        config.t0          = 0;
        config.tmax        = 200;
        config.dt          = 0.1;
        config.lambdas     = {0};
        config.gamma       = 1. / 7.;
        config.nu          = 1 / 20.;
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
        config.tmax                   = 3 * 365; // From 1st August 2016 to 1st August 2019
        config.dt                     = 0.1;
        config.lambdas                = {0.000002}; // TODO
        config.gamma                  = 1. / 7.;
        config.nu                     = 1 / 20.;
        config.I0s                    = {{0, 1000}}; // TODO
        config.R0s                    = {{0, 20000}}; // TODO
        config.total_populations      = {100000};
        config.first_season_start_day = -153; // First season starts at 1st August of the previous year
        config.season_peaks           = {5, 20, 28}; // Peak transmission rate at day January 24th of each season
        config.seasonality_rhos       = {0.8, 0.9, 0.7}; // Variation factor for each season
        config.seasonality_sigmas     = {40, 40, 60}; // Standard deviation for the Gaussian function for each season
        return config;
        break;
    }
    return config;
}

} // namespace Config
