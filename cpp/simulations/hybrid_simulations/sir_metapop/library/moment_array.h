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

#ifndef MOMENT_ARRAY_H
#define MOMENT_ARRAY_H

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>
#include "memilio/math/eigen.h"

/**
 * @brief A class template for a moment array.
 * @tparam NumInfectionStates Number of infection states in the model.
 * @tparam NumRegions Number of regions in the model; there is one index per infection state and region.
 * @tparam MaxIndex Maximum value each index can take; determines the maximum order of moments that can be stored.
 */
template <size_t NumInfectionStates, size_t NumRegions, size_t MaxIndex>
class MomentArray
{
public:
    using ArrayType = Eigen::Array<double, Eigen::Dynamic, 1>;

    MomentArray()
        : m_moments(ArrayType::Constant(std::pow(MaxIndex + 1, NumInfectionStates * NumRegions), 1, 0.))
    {
    }

    /**
     * @brief Takes an array of indices (one for each subpopulation/species i.e. one per region and infection state) and calculates the flat index of the corresponding moment.
     */
    size_t flatten_index(const std::array<int, NumInfectionStates * NumRegions>& indices) const
    {
        size_t index = 0;
        for (size_t i = 0; i < indices.size(); ++i) {
            index = index * (MaxIndex + 1) + indices[i];
        }
        return index;
    }

    /**
     * @brief Inversion function to flatten_index. Takes a flat index and calculates the multiindex of the corresponding moment.
     */
    std::array<int, NumInfectionStates * NumRegions> unflatten_index(size_t flat_index)
    {
        assert(flat_index < static_cast<size_t>(m_moments.rows()) &&
               "Flat index is bigger than number of moments for given order.");
        std::array<int, NumInfectionStates * NumRegions> indices{};
        for (size_t i = 0; i < indices.size(); ++i) {
            indices[NumInfectionStates * NumRegions - 1 - i] = flat_index % (MaxIndex + 1);
            flat_index /= (MaxIndex + 1);
        }
        return indices;
    }

    /**
     * @brief Returns reference to internally stored flat array.
     */
    ArrayType& moments()
    {
        return m_moments;
    }

    /**
     * @brief Returns a vector of the moments up to a given order.
     * @param[in] order Maximum order of moments to be returned.
     */
    std::vector<double> moments_up_to_order(int order)
    {
        std::vector<double> moments;
        for (size_t i = 0; i < static_cast<size_t>(m_moments.rows()); ++i) {
            auto multi_idx = unflatten_index(i);
            int sum        = 0;
            for (auto&& idx : multi_idx) {
                sum += idx;
            }
            if (sum <= order) {
                moments.push_back(m_moments[i]);
            }
        }
        return moments;
    }

    /**
     * @brief Returns a vector of the moment names up to a given order. Corresponds to moments calculated by moments_up_to_order.
     * @param[in] order Maximum order of moment names to be returned.
     */
    std::vector<std::string> names_up_to_order(int order)
    {
        std::vector<std::string> names;
        for (size_t i = 0; i < static_cast<size_t>(m_moments.rows()); ++i) {
            auto multi_idx = unflatten_index(i);
            int sum        = 0;
            for (auto&& idx : multi_idx) {
                sum += idx;
            }
            if (sum <= order) {
                std::string moment_name = "M";
                for (auto&& idx : multi_idx) {
                    moment_name += std::to_string(idx);
                }
                names.push_back(moment_name);
            }
        }
        return names;
    }

    /**
     * @brief Returns the array entry (moment value) given a multiindex.
     */
    double& operator[](const std::array<int, NumInfectionStates * NumRegions>& indices)
    {
        return m_moments[flatten_index(indices)];
    }

    /**
     * @brief Returns a vector of the moment names in the same order as they are stored in m_moments.
     */
    std::vector<std::string> get_names()
    {

        std::vector<std::string> names(m_moments.rows());
        for (size_t index = 0; index < static_cast<size_t>(m_moments.rows()); ++index) {
            std::string moment_name = "M";
            auto multi_idx          = unflatten_index(index);
            for (auto&& i : multi_idx) {
                moment_name += std::to_string(i);
            }
            names[index] = moment_name;
        }
        return names;
    }

private:
    ArrayType m_moments;
};

#endif //MOMENT_ARRAY_H
