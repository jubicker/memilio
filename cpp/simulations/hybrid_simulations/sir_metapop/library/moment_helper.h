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
#ifndef MOMENT_HELPER_H
#define MOMENT_HELPER_H

#include "simulations/hybrid_simulations/sir_metapop/config/config.h"
#include "moment_array.h"
#include "smm_moments/model.h"
#include "smm_moments/parameters.h"
#include "ode_sir/infection_state.h"
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <cstddef>
#include <numeric>

namespace moment_helper
{

/**
 * @brief Initialize moment model for given config.
 * @tparam NumRegions Number of regions.
 * @tparam ClosureOrder Order used for zero cumulant closure.
 * @param[in] expected_values_init Initial expected values.
 * @param[in] moment_init Initial moment values.
 * @param[in] config Configuration
 */
template <size_t NumRegions, size_t ClosureOrder>
mio::smm_moments::Model<NumRegions, ClosureOrder>
initialize_model(Eigen::Array<double, Eigen::Dynamic, 1>& expected_values_init,
                 Eigen::Array<double, Eigen::Dynamic, 1>& moments_init, const Config::Config& config,
                 typename mio::smm_moments::Model<NumRegions, ClosureOrder>::ClosureFunctionType closure_func)
{
    mio::smm_moments::Model<NumRegions, ClosureOrder> model(closure_func);
    // Check whether initial expected values and moments have the correct size
    assert(expected_values_init.rows() == NumRegions * static_cast<size_t>(mio::osir::InfectionState::Count) &&
           "Initial expected values do not have correct size");
    assert(moments_init.rows() == model.moments.moments().rows() && "Initial moments do not have correct size");

    // Set spatial transition rates
    for (auto& rate : config.transition_rates) {
        model.parameters.template get<mio::smm_moments::TransitionRate>()[{rate.status, rate.from, rate.to}] =
            rate.factor;
    }

    // Set recovery rate
    model.parameters.template get<mio::smm_moments::RecoveryRate>() = config.gamma;

    // Set immunity loss rate if > 0
    if (config.nu > 0) {
        model.parameters.template get<mio::smm_moments::ImmunityLossRate>() = config.nu;
    }

    // Set seasonality parameters if there are any
    if (config.season_peaks.size() > 0) {
        for (size_t season = 0; season < config.seasonality_rhos.size(); ++season) {
            model.parameters.template get<mio::smm_moments::SeasonalityRho>().push_back(
                config.seasonality_rhos[season]);
            model.parameters.template get<mio::smm_moments::SeasonalitySigma>().push_back(
                config.seasonality_sigmas[season]);
            model.parameters.template get<mio::smm_moments::SeasonalityPeak>().push_back(config.season_peaks[season]);
        }
        model.parameters.template get<mio::smm_moments::StartDay>()            = config.first_season_start_day;
        model.parameters.template get<mio::smm_moments::FirstSeasonStartDay>() = config.first_season_start_day;
    }

    for (size_t r = 0; r < NumRegions; ++r) {
        // Set transmission rates
        model.parameters.template get<mio::smm_moments::TransmissionRate>()[mio::regions::Region(r)] =
            config.lambdas[r];

        // Set initial expected values
        model.populations[{mio::regions::Region(r), mio::osir::InfectionState::Susceptible}] =
            expected_values_init[r * static_cast<size_t>(mio::osir::InfectionState::Count) +
                                 static_cast<size_t>(mio::osir::InfectionState::Susceptible)];
        model.populations[{mio::regions::Region(r), mio::osir::InfectionState::Infected}] =
            expected_values_init[r * static_cast<size_t>(mio::osir::InfectionState::Count) +
                                 static_cast<size_t>(mio::osir::InfectionState::Infected)];
        model.populations[{mio::regions::Region(r), mio::osir::InfectionState::Recovered}] =
            expected_values_init[r * static_cast<size_t>(mio::osir::InfectionState::Count) +
                                 static_cast<size_t>(mio::osir::InfectionState::Recovered)];
    }
    // Set initial moments
    model.moments.moments() = moments_init;
    return model;
}

/**
 * @brief Split line by comma.
 */
void split_line(std::string string, std::vector<double>* row)
{
    std::vector<std::string> strings;
    boost::split(strings, string, boost::is_any_of(","));
    std::transform(strings.begin(), strings.end(), std::back_inserter(*row), [&](std::string s) {
        return std::stod(s);
    });
}

/**
 * @brief Convert Moment string to the corresponding multi-index.
 */
template <size_t NumInfectionStates, size_t NumRegions>
std::array<int, NumInfectionStates * NumRegions> moment_to_indices(std::string moment)
{
    std::array<int, NumInfectionStates * NumRegions> moment_index{};
    for (size_t i = 0; i < moment_index.size(); ++i) {
        moment_index[i] = moment[i + 1] - '0'; // convert char to int
    }
    return moment_index;
}

/**
 * @brief Read expected values from csv timeseries file.
 * @param[in] file_path Csv file.
 * @param[in] time Time point from which values should be taken.
 * @return Array with expected values.
 */
Eigen::Array<double, Eigen::Dynamic, 1> read_expected_values(const std::string& file_path, double time)
{
    const boost::filesystem::path f = file_path;
    if (!boost::filesystem::exists(f)) {
        mio::log_error("Cannot read in data. Expected values file does not exist.");
    }
    // File pointer
    std::fstream fin_f;
    // Open file
    fin_f.open(file_path, std::ios::in);
    std::string line_f;
    // Read the titles from file
    std::getline(fin_f, line_f);
    std::vector<std::string> titles;
    boost::split(titles, line_f, boost::is_any_of(","));
    uint32_t col_count = titles.size();
    // First column is time, so we subtract 1
    Eigen::Array<double, Eigen::Dynamic, 1> expected_values_init(col_count - 1);
    // Read until the correct time is found
    std::vector<double> row;
    while (std::getline(fin_f, line_f)) {
        row.clear();
        // read columns in this row
        split_line(line_f, &row);
        if (row[0] == time) {
            for (size_t i = 1; i < col_count; ++i) {
                expected_values_init[i - 1] = row[i];
            }
            break;
        }
    }
    fin_f.close();
    return expected_values_init;
}

/**
 * @brief Read moments from csv timeseries file.
* @param[in] file_path Csv file.
 * @param[in] time Time point from which values should be taken.
 * @return Array with moments.
 */
template <size_t NumRegions, size_t ClosureOrder>
Eigen::Array<double, Eigen::Dynamic, 1> read_moments(const std::string& file_path, double time)
{
    const boost::filesystem::path f = file_path;
    if (!boost::filesystem::exists(f)) {
        mio::log_error("Cannot read in data. Moment file does not exist.");
    }
    // File pointer
    std::fstream fin_f;
    // Open file
    fin_f.open(file_path, std::ios::in);
    std::string line_f;
    // Read the titles from file
    std::getline(fin_f, line_f);
    //line_hosp.erase(std::remove(line_hosp.begin(), line_hosp.end(), '\r'), line_hosp.end());
    std::vector<std::string> titles;
    boost::split(titles, line_f, boost::is_any_of(","));
    uint32_t col_count = 0;
    // Get title by index
    std::map<std::string, uint32_t> index = {};
    for (auto const& title : titles) {
        index.insert({title, col_count});
        col_count++;
    }
    MomentArray<static_cast<size_t>(mio::osir::InfectionState::Count), NumRegions, ClosureOrder> moments_array;
    // Read until the correct time is found
    std::vector<double> row;
    while (std::getline(fin_f, line_f)) {
        row.clear();
        // read columns in this row
        split_line(line_f, &row);
        if (row[index["Time"]] == time) {
            for (const auto& [key, _] : index) {
                if (key == "Time") {
                    continue;
                }
                // Get multi-index from string
                const auto moment_index =
                    moment_to_indices<static_cast<size_t>(mio::osir::InfectionState::Count), NumRegions>(key);
                auto order = std::accumulate(moment_index.begin(), moment_index.end(), 0);
                if (static_cast<size_t>(order) >= ClosureOrder) {
                    continue;
                }
                // Get flat index from multi-index
                size_t flat_index                   = moments_array.flatten_index(moment_index);
                moments_array.moments()[flat_index] = row[index[key]];
            }
            break;
        }
    }
    fin_f.close();
    return moments_array.moments();
}

} // namespace moment_helper
#endif // MOMENT_HELPER_H
