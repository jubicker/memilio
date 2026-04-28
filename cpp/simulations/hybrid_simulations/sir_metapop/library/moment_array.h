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

#include <Eigen/src/Core/Array.h>
#include <boost/math/tools/mp.hpp>
#include <cmath>
#include <cstddef>
#include <map>
#include <numeric>
#include <string>
#include <vector>
#include "memilio/math/eigen.h"

/**
 * @brief A class template for a moment array.
 * @tparam NumInfectionStates Number of infection states in the model.
 * @tparam NumRegions Number of regions in the model; there is one index per infection state and region.
 * @tparam MaxOrder Maximum order of moments that can be stored.
 */
template <size_t NumInfectionStates, size_t NumRegions, size_t MaxOrder>
class MomentArray
{
public:
    using ArrayType  = Eigen::Array<double, Eigen::Dynamic, 1>;
    using MultiIndex = std::array<int, NumInfectionStates * NumRegions>;

    MomentArray()
    {
        // Once add all MultiIndices with order <= MaxOrder to m_indices
        MultiIndex idx{};
        fill_indices(0, 0, idx);
        // Initialize moments with correct size
        m_moments = ArrayType::Zero(m_indices.size());
        // Hash lookup for flat index
        for (size_t i = 0; i < m_indices.size(); ++i) {
            m_lookup[m_indices[i]] = i;
        }
        // Calculate orders
        m_orders.reserve(m_indices.size());
        for (auto& m : m_indices) {
            m_orders.push_back(std::accumulate(m.begin(), m.end(), 0));
        }
        // Fill moment names
        m_names.reserve(m_indices.size());
        fill_moment_names();
    }

    /**
     * @brief Takes an array of indices (One per region and infection state) and returns the flat index of the corresponding moment.
     */
    size_t flatten_index(const MultiIndex& idx) const
    {
        return m_lookup.at(idx);
    }

    /**
     * @brief Inversion function to flatten_index. Takes a flat index and returns the multi-index of the corresponding moment.
     */
    MultiIndex unflatten_index(size_t flat_index) const
    {
        assert(flat_index < static_cast<size_t>(m_moments.rows()) &&
               "Flat index is bigger than number of moments for given order.");
        return m_indices[flat_index];
    }

    /**
     * @brief Returns reference to internally stored flat array.
     */
    ArrayType moments() const
    {
        return m_moments;
    }

    ArrayType& moments()
    {
        return m_moments;
    }

    /**
     * @brief Returns MultiIndex vector for given moments.
     */
    const std::vector<MultiIndex>& get_indices() const
    {
        return m_indices;
    }

    /*
     * @brief Returns MultiIndex at position flat_index. 
     */
    MultiIndex& index(size_t flat_index) const
    {
        return m_indices[flat_index];
    }

    /**
     * @brief Returns MultiIndex strings for given moments.
     */
    const std::vector<std::string>& get_names() const
    {
        return m_names;
    }

    /**
     * @brief Returns orders for given moments.
     */
    const std::vector<size_t>& get_order() const
    {
        return m_orders;
    }

    /**
     * @brief Returns order of moment at position flat_index.
     */
    size_t order(size_t flat_index) const
    {
        return m_orders[flat_index];
    }

    /**
     * @brief Returns moment values.
     */
    ArrayType get_values() const
    {
        return m_moments;
    }

    /**
     * @brief Returns the array entry (moment value) given a multi-index.
     */
    double& operator[](const std::array<int, NumInfectionStates * NumRegions>& indices)
    {
        return m_moments[flatten_index(indices)];
    }

private:
    /**
     * @brief Recursively adds all MultiIndices with order <= MaxOrder to m_indices vector.
     * @param[in] pos Current position in MultiInddex that is to fill.
     * @param[in] sum Sum of the MultiIndex positions 0,...,pos-1 that have already been set.
     * @param[in, out] idx Current MultiIndex whose values 0,..,pos-1 have already been filled. Is added to m_indices when fully filled.
     */
    void fill_indices(size_t pos, size_t sum, MultiIndex& idx)
    {
        if (pos == NumInfectionStates * NumRegions) {
            // full multi-index
            m_indices.push_back(idx);
            return;
        }

        const size_t maxAllowed = MaxOrder - sum;
        for (size_t v = 0; v <= maxAllowed; ++v) {
            idx[pos] = v;
            fill_indices(pos + 1, sum + v, idx);
        }
    }

    /**
     * @brief Returns multi-indices as string for given moments.
     */
    void fill_moment_names()
    {
        for (const auto& m : m_indices) {
            std::string moment_name = "M";
            for (auto&& idx : m) {
                moment_name += std::to_string(idx);
            }
            m_names.push_back(moment_name);
        }
    }

    std::vector<MultiIndex> m_indices; ///< Multi-index of all moments with order <= MaxOrder.
    std::vector<std::string> m_names; ///< m_indices as string vector.
    std::vector<size_t> m_orders; ///< Orders of all moments considered.
    std::map<MultiIndex, size_t>
        m_lookup; ///< Map that maps a moment's multi-index to its corresponding flat index in m_moments.
    ArrayType m_moments; ///< Array containing all moment values.
};

#endif //MOMENT_ARRAY_H
