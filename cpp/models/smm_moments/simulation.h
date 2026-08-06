/* 
* Copyright (C) 2020-2025 German Aerospace Center (DLR-SC)
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
#ifndef MOMENTS_SIMULATION_H
#define MOMENTS_SIMULATION_H

#include "memilio/utils/compiler_diagnostics.h"
#include "ode_sir/infection_state.h"
#include "smm_moments/model.h"
#include "memilio/compartments/simulation_base.h"
#include "memilio/config.h"
#include "memilio/math/stepper_wrapper.h"
#include "memilio/utils/time_series.h"
#include <cassert>
#include <cstddef>
#include <limits>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

namespace mio
{

namespace smm_moments
{

/// @brief The default integrator used by Simulation.
template <typename FP>
using DefaultIntegratorCore = mio::ControlledStepperWrapper<FP, boost::numeric::odeint::runge_kutta_cash_karp54>;

template <size_t NumRegions, size_t ClosureOrder>
class Simulation
    : public mio::details::SimulationBase<ScalarType, Model<NumRegions, ClosureOrder>, mio::OdeIntegrator<ScalarType>>
{
public:
    using Base =
        mio::details::SimulationBase<ScalarType, Model<NumRegions, ClosureOrder>, mio::OdeIntegrator<ScalarType>>;

    /**
     * @brief Setup the simulation with an ODE solver.
     * @param[in] model An instance of the moment model
     * @param[in] t0 Start time.
     * @param[in] dt Initial step size of integration
     */
    Simulation(Model<NumRegions, ClosureOrder> const& model, ScalarType t0 = 0., ScalarType dt = 0.1)
        : Base(model, std::make_unique<DefaultIntegratorCore<ScalarType>>(), t0, dt)
    {
    }

    /**
     * @brief Run the simulation up to a given time.
     * The time tmax must be greater than `get_result().get_last_time_point()`, which is used as the starting point. The initial value is `get_result().get_last_value()`.
     * @param[in] tmax Next stopping point of the simulation.
     * @return The simulation result at tmax.
     */
    Eigen::Ref<Eigen::VectorX<ScalarType>> advance(ScalarType tmax)
    {
        return Base::advance(
            [this](auto&& y, auto&& t, auto&& dydt) {
                Base::get_model().eval_right_hand_side(y, y, t, dydt);
            },
            tmax, Base::get_result());
    }

    /**
     * @brief Extracts time series of only the expected values from simulation result.
     * @return Time series of expected values.
     */
    mio::TimeSeries<double> get_expected_values_time_series() const
    {
        mio::TimeSeries<double> expected_values_ts(static_cast<size_t>(mio::osir::InfectionState::Count) * NumRegions);
        for (auto t = 0; t < Base::get_result().get_num_time_points(); ++t) {
            auto y = Base::get_result().get_value(t);
            expected_values_ts.add_time_point(Base::get_result().get_time(t),
                                              y.head(Base::get_model().populations.get_num_compartments()));
        }
        return expected_values_ts;
    }

    /**
     * @brief Extracts time series of all moments up to given order from simulation result.
     * @param order Maximum order of moments to extract.
     * @return Pair of time series of moments and corresponding moment names.
     */
    std::pair<mio::TimeSeries<double>, std::vector<std::string>> get_moment_time_series(size_t order)
    {
        assert(order <= ClosureOrder);
        mio::TimeSeries<double> moment_ts(Base::get_model().moments.get_indices().size());
        // Add first moments and fill name vector
        size_t index               = 0;
        Eigen::VectorXd moment_vec = Eigen::VectorXd::Zero(moment_ts.get_num_elements());
        // Add remaining time points
        for (auto t = 0; t < Base::get_result().get_num_time_points(); ++t) {
            index      = 0;
            moment_vec = Eigen::VectorXd::Zero(moment_ts.get_num_elements());
            auto y     = Base::get_result().get_value(t);
            for (size_t i = 0; i < y.size() - Base::get_model().populations.get_num_compartments(); i++) {
                auto multi_idx = Base::get_model().moments.unflatten_index(i);
                int sum        = std::accumulate(multi_idx.begin(), multi_idx.end(), 0);
                if (sum <= int(order)) {
                    moment_vec(index) = y[Base::get_model().populations.get_num_compartments() + i];
                    index++;
                    if (index == static_cast<size_t>(moment_vec.size())) {
                        break;
                    }
                }
            }
            moment_ts.add_time_point(Base::get_result().get_time(t), moment_vec);
        }
        return std::make_pair(moment_ts, Base::get_model().moments.get_names());
    }

    /**
     * @brief Extracts the last value of all variances from simulation result.
     * @return Vector of all variances.
     */
    std::vector<ScalarType> get_last_vars()
    {
        std::vector<ScalarType> vars(static_cast<size_t>(osir::InfectionState::Count) * NumRegions);
        // Add last variances
        auto y = Base::get_result().get_last_value();
        for (size_t i = 0; i < y.size() - Base::get_model().populations.get_num_compartments(); i++) {
            auto multi_idx = Base::get_model().moments.unflatten_index(i);
            bool is_var    = std::count(multi_idx.begin(), multi_idx.end(), 2) == 1 &&
                          std::count(multi_idx.begin(), multi_idx.end(), 0) == multi_idx.size() - 1;
            if (!is_var) {
                continue;
            }
            size_t index = std::distance(multi_idx.begin(), std::find(multi_idx.begin(), multi_idx.end(), 2));
            vars[index]  = y[Base::get_model().populations.get_num_compartments() + i];
        }

        return vars;
    }

    /**
     * @brief Extracts the last gradient of all variances from simulation result.
     * @return Vector of all variance gradients.
     */
    std::vector<ScalarType> get_last_var_gradients()
    {
        std::vector<ScalarType> vars_gradient(static_cast<size_t>(osir::InfectionState::Count) * NumRegions);
        const bool has_gradient = Base::get_result().get_num_time_points() > 1;
        auto last_tp            = Base::get_result().get_last_time();
        auto y_last             = Base::get_result().get_last_value();
        auto second_last_tp_index = has_gradient ? Base::get_result().get_num_time_points() - 2 : 0;
        auto second_last_tp       = has_gradient ? Base::get_result().get_time(second_last_tp_index) : last_tp;
        auto y_second_last        = has_gradient ? Base::get_result().get_value(second_last_tp_index) : y_last;
        for (size_t i = 0; i < y_last.size() - Base::get_model().populations.get_num_compartments(); i++) {
            auto multi_idx = Base::get_model().moments.unflatten_index(i);
            bool is_var    = std::count(multi_idx.begin(), multi_idx.end(), 2) == 1 &&
                          std::count(multi_idx.begin(), multi_idx.end(), 0) == multi_idx.size() - 1;
            if (!is_var) {
                continue;
            }
            size_t index = std::distance(multi_idx.begin(), std::find(multi_idx.begin(), multi_idx.end(), 2));
            // If there is only one time point, the gradient is set to max double value to ensure that it is above any reasonable threshold for switching conditions
            vars_gradient[index] = has_gradient ? (y_last[Base::get_model().populations.get_num_compartments() + i] -
                                                   y_second_last[Base::get_model().populations.get_num_compartments() + i]) /
                                                      (last_tp - second_last_tp)
                                                 : std::numeric_limits<double>::max();
        }

        return vars_gradient;
    }

    /**
     * @brief Extracts the last gradient of all means from simulation result.
     * @return Vector of all mean gradients.
     */
    std::vector<ScalarType> get_last_mean_gradients()
    {
        std::vector<ScalarType> means_gradient(static_cast<size_t>(osir::InfectionState::Count) * NumRegions);
        const bool has_gradient = Base::get_result().get_num_time_points() > 1;
        auto last_tp            = Base::get_result().get_last_time();
        auto y_last             = Base::get_result().get_last_value();
        auto second_last_tp_index = has_gradient ? Base::get_result().get_num_time_points() - 2 : 0;
        auto second_last_tp       = has_gradient ? Base::get_result().get_time(second_last_tp_index) : last_tp;
        auto y_second_last        = has_gradient ? Base::get_result().get_value(second_last_tp_index) : y_last;
        for (size_t i = 0; i < Base::get_model().populations.get_num_compartments(); i++) {
            // If there is only one time point, the gradient is set to max double value to ensure that it is above any reasonable threshold for switching conditions
            means_gradient[i] = has_gradient ? (y_last[i] - y_second_last[i]) / (last_tp - second_last_tp)
                                              : std::numeric_limits<double>::max();
        }

        return means_gradient;
    }

    /**
     * @brief Extracts the last value of all means from simulation result.
     * @return Vector of all means.
     */
    std::vector<ScalarType> get_last_means()
    {
        auto y    = Base::get_result().get_last_value();
        auto head = y.head(Base::get_model().populations.get_num_compartments());
        std::vector<ScalarType> means(head.data(), head.data() + head.size());

        return means;
    }
};

} // namespace smm_moments

} // namespace mio

#endif // MOMENTS_SIMULATION_H
