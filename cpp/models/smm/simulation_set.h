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

#ifndef MIO_SMM_SIMULATION_SET_H
#define MIO_SMM_SIMULATION_SET_H

#include "memilio/config.h"
#include "memilio/timer/basic_timer.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/random_number_generator.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/mioomp.h"
#include "ode_sir/infection_state.h"
#include "smm/simulation.h"
#include "simulations/hybrid_simulations/sir_metapop/library/moment_array.h"
#include "memilio/data/analyze_result.h"
#include <cstddef>
#include <numeric>
#include <sys/types.h>
#include <vector>

namespace mio
{
namespace smm
{

/**
 * @brief A class performing multiple simulations of mio::smm::Model.
 * @tparam regions The number of regions.
 * @tparam Status The used infection state enum
 * @tparam MaxMomentOrder The maximum order up to which moments are calculated as result.
 */
template <size_t regions, class Status, size_t MaxMomentOrder>
class SimulationSet
{
public:
    using Simulation = smm::Simulation<ScalarType, regions, Status>;
    using Model      = Simulation::Model;

    /**
     * @brief Set up multiple simulations for SMM.
     * @param[in] num_runs Number of simulations.
     * @param[in] model Model which is simulated multiple times.
     * @param[in] t0 Start time.
     * @param[in] dt Interpolation step size and initial simulation time step.
     */
    SimulationSet(size_t num_runs, const Model& model, ScalarType t0, ScalarType dt)
        : m_models(num_runs, model)
        , m_means(static_cast<size_t>(Status::Count) * regions)
        , m_moments(1)
        , m_results(num_runs, TimeSeries<ScalarType>(static_cast<size_t>(Status::Count) * regions))
        , m_dt(dt)
        , m_t(t0)
        , m_t_index(0)
        , m_sim_time(num_runs, 0.)
        , last_results(Eigen::Matrix<double, Eigen::Dynamic, static_cast<size_t>(Status::Count) * regions>::Zero(
              num_runs, static_cast<size_t>(Status::Count) * regions))
        , means(Eigen::Matrix<double, 1, static_cast<size_t>(Status::Count) * regions>::Zero(
              1, static_cast<size_t>(Status::Count) * regions))
        , centered_values(Eigen::Matrix<double, Eigen::Dynamic, static_cast<size_t>(Status::Count) * regions>::Zero(
              num_runs, static_cast<size_t>(Status::Count) * regions))
    {
        for (auto& m : m_models) {
            m.get_rng() = mio::RandomNumberGenerator();
            m_sims.push_back(Simulation(m, t0, dt));
        }
        m_moments = TimeSeries<double>(m_mom_array.get_indices().size());

        // Add initial values to results
        auto& sim_ts = m_sims[0].get_result();
#ifdef MEMILIO_ENABLE_OPENMP
#pragma omp parallel for
#endif
        for (size_t run = 0; run < m_sims.size(); ++run) {
            if (m_results[run].get_num_time_points() == 0) {
                m_results[run].add_time_point(sim_ts.get_time(0), sim_ts.get_value(0));
            }
        }
        update_means_and_moments();
    }

    /**
     * @brief Advance all simulations until tmax.
     * Additionally, the advance time per run is saved and total time is increased by runtime of this function.
     * @param[in] tmax Simulation end time point.
     */
    void advance(double tmax)
    {

// Run simulations
#ifdef MEMILIO_ENABLE_OPENMP
#pragma omp parallel for
#endif
        for (size_t run = 0; run < m_sims.size(); ++run) {
            timing::BasicTimer timer;
            timer.start();
            m_sims[run].advance(tmax);
            timer.stop();
            m_sim_time[run] += mio::timing::time_in_seconds(timer.get_elapsed_time());
        }

        m_t = tmax;

        calculate_outputs();
    }

    /**
     * @brief Calculate interpolated results, mean and moment time series.
     */
    void calculate_outputs()
    {
        double t = 0;
        if (m_means.get_num_time_points() > 0) {
            t = m_means.get_last_time();
        }
        // Time steps for interpolation
        int num_steps = static_cast<int>((m_t - t) / m_dt);
        if (num_steps > 0) {

// Interpolate results
#ifdef MEMILIO_ENABLE_OPENMP
#pragma omp parallel for
#endif
            for (size_t run = 0; run < m_sims.size(); ++run) {
                auto& sim_ts         = m_sims[run].get_result();
                auto interpolated_ts = sim_ts;
                // Remove individual simulation result time series
                sim_ts = mio::TimeSeries<double>(interpolated_ts.get_num_elements());
                while (interpolated_ts.get_num_time_points() > 1) {
                    if (m_results[run].get_num_time_points() == 0 ||
                        m_results[run].get_last_time() < interpolated_ts.get_time(0)) {
                        m_results[run].add_time_point(interpolated_ts.get_time(0), interpolated_ts.get_value(0));
                    }
                    interpolated_ts.remove_time_point(0);
                }
                m_results[run].add_time_point(interpolated_ts.get_last_time(), interpolated_ts.get_last_value());
                sim_ts.add_time_point(interpolated_ts.get_last_time(), interpolated_ts.get_last_value());
            }

            // Fill moment and means time series
            update_means_and_moments();
        }
    }

    /**
     * @brief Get all simulations.
     */
    std::vector<Simulation>& get_simulations()
    {
        return m_sims;
    }
    const std::vector<Simulation>& get_simulations() const
    {
        return m_sims;
    }

    /**
     * @brief Get interpolated time series of all runs.
     */
    std::vector<TimeSeries<double>>& get_result()
    {
        return m_results;
    }
    const std::vector<TimeSeries<double>>& get_result() const
    {
        return m_results;
    }

    /**
     * @brief Get mean time series.
     */
    TimeSeries<double>& get_mean()
    {
        return m_means;
    }
    const TimeSeries<double>& get_mean() const
    {
        return m_means;
    }

    /**
     * @brief Get last mean value.
     */
    std::vector<ScalarType> get_last_means()
    {
        auto last_value = m_means.get_last_value();
        std::vector<ScalarType> last_means(last_value.data(), last_value.data() + last_value.size());
        return last_means;
    }

    /**
     * @brief Get moment time series.
     */
    TimeSeries<double>& get_moments()
    {
        return m_moments;
    }
    const TimeSeries<double>& get_moments() const
    {
        return m_moments;
    }

    /**
     * @brief Get last value of all variances.
     */
    std::vector<double> get_last_vars()
    {
        std::vector<double> vars(static_cast<size_t>(osir::InfectionState::Count) * regions);
        auto& names = m_mom_array.get_names();
#ifdef MEMILIO_ENABLE_OPENMP
#pragma omp parallel for
#endif
        for (size_t i = 0; i < names.size(); ++i) {
            bool is_var = std::count(names[i].begin(), names[i].end(), '2') == 1 &&
                          std::count(names[i].begin(), names[i].end(), '0') ==
                              static_cast<size_t>(osir::InfectionState::Count) * regions - 1;
            if (!is_var) {
                continue;
            }
            size_t index = std::distance(names[i].begin(), std::find(names[i].begin(), names[i].end(), '2')) - 1;
            vars[index]  = m_moments.get_last_value()[i];
        }
        return vars;
    }

    /**
     * @brief Get last gradient of all variances.
     */
    std::vector<double> get_last_var_gradients()
    {
        std::vector<double> vars_gradient(static_cast<size_t>(osir::InfectionState::Count) * regions);
        auto& names = m_mom_array.get_names();
#ifdef MEMILIO_ENABLE_OPENMP
#pragma omp parallel for
#endif
        for (size_t i = 0; i < names.size(); ++i) {
            bool is_var = std::count(names[i].begin(), names[i].end(), '2') == 1 &&
                          std::count(names[i].begin(), names[i].end(), '0') ==
                              static_cast<size_t>(osir::InfectionState::Count) * regions - 1;
            if (!is_var) {
                continue;
            }
            size_t index = std::distance(names[i].begin(), std::find(names[i].begin(), names[i].end(), '2')) - 1;
            auto last_tp = m_moments.get_last_time();
            if (m_moments.get_num_time_points() > 1) {
                auto second_last_tp_index = m_moments.get_num_time_points() - 2;
                auto second_last_tp       = m_moments.get_time(second_last_tp_index);
                auto y_second_last        = m_moments.get_value(second_last_tp_index);
                auto y_last               = m_moments.get_last_value();
                vars_gradient[index]      = (y_last[i] - y_second_last[i]) / (last_tp - second_last_tp);
            }
            else {
                vars_gradient[index] = 0.;
            }
        }
        return vars_gradient;
    }

    /**
     * @brief Get moment names.
     */
    std::vector<std::string> get_moment_names()
    {
        return m_mom_array.get_names();
    }
    const std::vector<std::string> get_moment_names() const
    {
        return m_mom_array.get_names();
    }

    const std::vector<std::array<int, static_cast<size_t>(Status::Count) * regions>>& get_moment_indices() const
    {
        return m_mom_array.get_indices();
    }

    const std::vector<size_t>& get_moment_orders() const
    {
        return m_mom_array.get_order();
    }

    /**
     * @brief Get run times of all runs.
     */
    std::vector<double>& get_sim_times()
    {
        return m_sim_time;
    }
    const std::vector<double>& get_sim_times() const
    {
        return m_sim_time;
    }

    /**
     * @brief Get maximum order of moments.
     */
    size_t get_order() const
    {
        return MaxMomentOrder;
    }

    /**
     * @brief Advance current time to t.
     * @param[in] t New current time.
     */
    void advance_t(double t)
    {
        if (t < m_t) {
            log_error("Current t of simulation set is set back. m_t is {} and is set to {}", m_t, t);
        }
        m_t = t;
    }

    /**
     * @brief Advance current time index by one.
     */
    void advance_t_index()
    {
        m_t_index += 1;
    }

    /**
     * @brief Calculate last mean and moment time point from simulation results.
     */
    void recalculate_last_means_and_moments()
    {
        means.setZero();
        last_results.setZero();
        centered_values.setZero();

#ifdef MEMILIO_ENABLE_OPENMP
#pragma omp parallel for
#endif
        // Copy simulation values to matrix
        for (size_t run = 0; run < m_results.size(); run++) {
            for (size_t r = 0; r < regions; ++r) {
                for (size_t s = 0; s < static_cast<size_t>(Status::Count); ++s) {
                    last_results(run, static_cast<size_t>(Status::Count) * r + s) =
                        m_results[run]
                            .get_last_value()[r * static_cast<size_t>(Status::Count) + static_cast<size_t>(Status(s))];
                }
            }
        }
        // Calculate means
        means = last_results.colwise().mean();
        if (m_means.get_num_time_points() == 0 || m_means.get_last_time() < m_results[0].get_last_time()) {
            m_means.add_time_point(m_results[0].get_last_time(), means);
        }
        else {
            m_means.get_last_value() = means;
        }

        // Calculate centered values used for calculation of central moments
        centered_values = last_results.rowwise() - means;

        auto& indices = m_mom_array.get_indices();

#ifdef MEMILIO_ENABLE_OPENMP
#pragma omp parallel for
#endif
        // Calculate moments
        for (auto& idx : indices) {
            m_mom_array[idx] = calculate_moment(centered_values, idx);
        }

        // Get only moments up to the given order
        auto moment_values   = m_mom_array.get_values();
        Eigen::VectorXd data = Eigen::VectorXd::Map(moment_values.data(), moment_values.size());
        if (m_moments.get_num_time_points() == 0 || m_moments.get_last_time() < m_results[0].get_last_time()) {
            m_moments.add_time_point(m_results[0].get_last_time(), data);
        }
        else {
            m_moments.get_last_value() = data;
        }
    }

    /**
     * @brief Set last mean values for given region for all infection states to 0.
     * @param[in] region Region whose mean values are set to 0.
     */
    void remove_means_for_region(size_t region)
    {
        if (region >= regions) {
            mio::log_error("Region doesn't exist. Region index is {} and number of modelled regions is {}", region,
                           regions);
            return;
        }
        for (size_t comp = 0; comp < static_cast<size_t>(Status::Count); ++comp) {
            m_means.get_last_value()[static_cast<size_t>(Status::Count) * region + comp] = 0;
        }
    }

    /**
     * @brief Set last moment values for moments that are fully or partially in given region to 0.
     * @param[in] region Region whose moments values are set to 0.
     */
    void remove_moments_for_region(size_t region)
    {
        auto& indices = m_mom_array.get_indices();

#ifdef MEMILIO_ENABLE_OPENMP
#pragma omp parallel for
#endif
        for (size_t i = 0; i < indices.size(); ++i) {
            size_t order_region = std::accumulate(indices[i].begin() + static_cast<size_t>(Status::Count) * region,
                                                  indices[i].begin() + static_cast<size_t>(Status::Count) * region +
                                                      static_cast<size_t>(Status::Count),
                                                  0);
            if (order_region > 0) {
                m_moments.get_last_value()[i] = 0;
            }
        }
    }

private:
    /**
     * @brief Calculate mean and moment time series from simulation results.
     */
    void update_means_and_moments()
    {
        for (auto t = m_t_index; t < m_results[0].get_num_time_points(); ++t) {
            means.setZero();
            last_results.setZero();
            centered_values.setZero();
            if (m_means.get_num_time_points() >= t + 1 && m_means.get_time(t) == m_results[0].get_time(t)) {
                log_warning("Mean time series already has time point t={}.", m_means.get_time(t));
                continue;
            }
#ifdef MEMILIO_ENABLE_OPENMP
#pragma omp parallel for
#endif
            // Copy simulation values to matrix
            for (size_t run = 0; run < m_results.size(); run++) {
                for (size_t r = 0; r < regions; ++r) {
                    for (size_t s = 0; s < static_cast<size_t>(Status::Count); ++s) {
                        last_results(run, static_cast<size_t>(Status::Count) * r + s) = m_results[run].get_value(
                            t)[r * static_cast<size_t>(Status::Count) + static_cast<size_t>(Status(s))];
                    }
                }
            }
            // Calculate means
            means = last_results.colwise().mean();
            m_means.add_time_point(m_results[0].get_time(t), means);

            // Calculate centered values used for calculation of central moments
            centered_values = last_results.rowwise() - means;

            auto& indices = m_mom_array.get_indices();
#ifdef MEMILIO_ENABLE_OPENMP
#pragma omp parallel for
#endif
            // Calculate moments
            for (auto& idx : indices) {
                m_mom_array[idx] = calculate_moment(centered_values, idx);
            }
            // Get only moments up to the given order
            auto moment_values   = m_mom_array.get_values();
            Eigen::VectorXd data = Eigen::VectorXd::Map(moment_values.data(), moment_values.size());
            m_moments.add_time_point(m_results[0].get_time(t), data);
            m_t_index += 1;
        }
    }

    /**
    * @brief Calculates moment defined by indices for the given realizations of a stochastic variable.
    * @param[in] values The realizations of the stochastic variable (i.e. S, I, R for each region).
    * @param[in] indices The indices defining the moment to be calculated.
    */
    double calculate_moment(
        const Eigen::Matrix<ScalarType, Eigen::Dynamic, static_cast<size_t>(Status::Count) * regions>& centered_vals,
        const std::array<int, static_cast<size_t>(Status::Count) * regions>& indices)
    {
        double moment = 0.0;
        for (int i = 0; i < centered_vals.rows(); ++i) {
            double summand = 1.0;
            for (size_t r = 0; r < regions; ++r) {
                for (size_t s = 0; s < static_cast<size_t>(Status::Count); ++s) {
                    summand *= std::pow(centered_vals(i, r * static_cast<size_t>(Status::Count) + s),
                                        indices[r * static_cast<size_t>(Status::Count) + s]);
                }
            }
            moment += summand;
        }
        moment /= static_cast<double>(centered_vals.rows());
        return moment;
    }

    std::vector<Model> m_models; ///< Models used for simulation.
    std::vector<Simulation> m_sims; ///< SMM simulations.
    TimeSeries<double> m_means; ///< Time series of means.
    TimeSeries<double> m_moments; ///< Time series of all moments up to MaxMomentOrder.
    std::vector<TimeSeries<ScalarType>> m_results; ///< Interpolated simulation results.
    double m_dt; ///< Interpolation time step.
    double m_t; ///< Current time.
    int m_t_index; ///< Current result time point index.
    std::vector<double> m_sim_time; ///< Simulation time per run.
    MomentArray<static_cast<size_t>(Status::Count), regions, MaxMomentOrder>
        m_mom_array{}; ///< Moment array used to calculate moments and their names.

    Eigen::Matrix<ScalarType, Eigen::Dynamic, static_cast<size_t>(Status::Count) * regions>
        last_results; ///< Helper matrix used to calculate moments.
    Eigen::Matrix<ScalarType, 1, static_cast<size_t>(Status::Count) * regions>
        means; ///< Helper matrix used to calculate means and moments.
    Eigen::Matrix<ScalarType, Eigen::Dynamic, static_cast<size_t>(Status::Count) * regions>
        centered_values; ///< Helper matrix used to calculate moments.
};

} // namespace smm
} // namespace mio

#endif // MIO_SMM_SIMULATION_SET_H
