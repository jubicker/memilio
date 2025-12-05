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
#include "memilio/utils/time_series.h"
#include "smm/simulation.h"
#include "simulations/hybrid_simulations/sir_metapop/library/moment_array.h"
#include <cstddef>
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
    using Simulation = smm::Simulation<regions, Status>;
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
        , m_sim_time(num_runs, 0.)
        , m_advance_time(0.)
    {
        for (auto& m : m_models) {
            m.get_rng().seed(m.get_rng().generate_seeds());
            m_sims.push_back(Simulation(m, t0, dt));
        }
        m_moments      = TimeSeries<double>(m_mom_array.names_up_to_order(MaxMomentOrder).size());
        m_moment_names = m_mom_array.names_up_to_order(MaxMomentOrder);
    }

    /**
     * @brief Advance all simulations until tmax and calculate interpolated results, mean and moment time series.
     * Additionally, the advance time per run and the total time of the function are saved.
     * @param[in] tmax Simulation end time point.
     */
    void advance(double tmax)
    {
        timing::BasicTimer total_timer;
        total_timer.start();
        // Run simulations
#ifdef MEMILIO_ENABLE_OPENMP
// Parallel
#pragma omp parallel for
        for (size_t run = 0; run < m_sims.size(); ++run) {
            timing::BasicTimer timer;
            timer.start();
            m_sims[run].advance(tmax);
            timer.stop();
            m_sim_time[run] = timer.get_elapsed_time();
        }
#else
        // Serial
        for (size_t run = 0; run < m_sims.size(); ++run) {
            timing::BasicTimer timer;
            timer.start();
            m_sims[run].advance(tmax);
            timer.stop();
            m_sim_time[run] = timer.get_elapsed_time();
        }
#endif

        // Time steps for interpolation
        int num_steps = static_cast<int>(tmax / m_dt) + 1;
        std::vector<double> interpolation_tps(num_steps);
        for (int i = 0; i < num_steps; ++i) {
            interpolation_tps[i] = i * m_dt;
        }

// Interpolate results
#ifdef MEMILIO_ENABLE_OPENMP
// Parallel
#pragma omp parallel for
        for (size_t run = 0; run < m_sims.size(); ++run) {
            m_results[run] = interpolate_simulation_result(m_sims[run].get_result(), interpolation_tps);
        }
#else
        // Serial
        for (size_t run = 0; run < m_sims.size(); ++run) {
            m_results[run] = interpolate_simulation_result(m_sims[run].get_result(), interpolation_tps);
        }
#endif

        // Fill moment and means time series
        calculate_means();
        calculate_moments();

        // Stop total timer and return elapsed time
        total_timer.stop();
        m_advance_time += total_timer.get_elapsed_time();
    }

    /**
     * @brief Get time series of all runs.
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
     * @brief Get moment names.
     */
    std::vector<std::string> get_moment_names()
    {
        return m_moment_names;
    }
    const std::vector<std::string> get_moment_names() const
    {
        return m_moment_names;
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
     * @brief Get total advance time.
     */
    double get_advance_time() const
    {
        return m_advance_time;
    }

private:
    /**
     * @brief Calculate mean time series from simulation results.
     */
    void calculate_means()
    {
        const size_t num_elements = static_cast<size_t>(Status::Count) * regions;
        for (auto t = 0; t < m_results[0].get_num_time_points(); ++t) {
            Eigen::Matrix<ScalarType, num_elements, 1> means;
            means.setZero();
            if (m_means.get_num_time_points() >= t + 1 && m_means.get_time(t) == m_results[0].get_time(t)) {
                log_warning("Mean time series already has time point t={}.", m_means.get_time(t));
                continue;
            }
            // Sum up all values
            for (auto& res : m_results) {
                means += res.get_value(t);
            }
            // Average values by number of runs
            means /= m_results.size();
            m_means.add_time_point(m_results[0].get_time(t), means);
        }
    }

    /**
    * @brief Calculates moment defined by indices for the given realizations of a stochastic variable.
    * @param[in] values The realizations of the stochastic variable (i.e. S, I, R for each region).
    * @param[in] indices The indices defining the moment to be calculated.
    */
    double calculate_moment(
        const Eigen::Matrix<ScalarType, Eigen::Dynamic, static_cast<size_t>(Status::Count) * regions>& values,
        const std::array<int, static_cast<size_t>(Status::Count) * regions>& indices)
    {
        Eigen::Matrix<ScalarType, 1, static_cast<size_t>(Status::Count)* regions> means = values.colwise().mean();
        double moment                                                                   = 0.0;
        for (int i = 0; i < values.rows(); ++i) {
            double summand = 1.0;
            for (size_t r = 0; r < regions; ++r) {
                for (size_t s = 0; s < static_cast<size_t>(Status::Count); ++s) {
                    summand *= std::pow(values(i, r * static_cast<size_t>(Status::Count) + s) -
                                            means(r * static_cast<size_t>(Status::Count) + s),
                                        indices[r * static_cast<size_t>(Status::Count) + s]);
                }
            }
            moment += summand;
        }
        moment /= static_cast<double>(values.rows());
        return moment;
    }

    /**
     * @brief Calculate moment time series from simulation results.
     */
    void calculate_moments()
    {
        const size_t num_elements = static_cast<size_t>(Status::Count) * regions;
        for (int t = 0; t < m_results[0].get_num_time_points(); ++t) {
            if (m_moments.get_num_time_points() >= t + 1 && m_moments.get_time(t) == m_results[0].get_time(t)) {
                log_warning("Moment time series already has time point t={}.", m_moments.get_time(t));
                continue;
            }
            // Copy results to matrix
            Eigen::Matrix<double, Eigen::Dynamic, num_elements> result =
                Eigen::Matrix<double, Eigen::Dynamic, num_elements>::Zero(m_results.size(), num_elements);
            for (size_t run = 0; run < m_results.size(); run++) {
                for (size_t r = 0; r < regions; ++r) {
                    for (size_t s = 0; s < static_cast<size_t>(Status::Count); ++s) {
                        result(run, static_cast<size_t>(Status::Count) * r + s) = m_results[run].get_value(
                            t)[r * static_cast<size_t>(Status::Count) + static_cast<size_t>(Status(s))];
                    }
                }
            }

            std::array<int, num_elements> indices; // Vector with current indices
            std::function<void(int, int)> fill_moments =
                [&](int pos, int currentSum) { // pos: current position in indices, currentSum: sum of indices so far
                    if (pos == int(indices.size()) &&
                        std::accumulate(indices.begin(), indices.end(), 0) <=
                            int(MaxMomentOrder)) { // Position is at last index i.e. all indiced for the moment are filled
                        m_mom_array[indices] = calculate_moment(result, indices);
                        return;
                    }

                    int maxAllowedHere =
                        std::min(MaxMomentOrder, MaxMomentOrder - currentSum); //maximum allowed value for current index
                    for (int v = 0; v <= maxAllowedHere; ++v) { // Iterate over all values allowed for the current index
                        indices[pos] = v;
                        int newSum   = currentSum + v; // Increase sum by current index
                        fill_moments(
                            pos + 1,
                            newSum); // This triggers the next index to take all possible values given the value of the current index
                    }
                };

            fill_moments(0, 0); // Start with first index and sum 0
            // Get only moments up to the given order
            auto moment_values   = m_mom_array.moments_up_to_order(MaxMomentOrder);
            Eigen::VectorXd data = Eigen::VectorXd::Map(moment_values.data(), moment_values.size());
            m_moments.add_time_point(m_results[0].get_time(t), data);
        }
    }

    std::vector<Model> m_models; ///< Models used for simulation.
    std::vector<Simulation> m_sims; ///< SMM simulations.
    TimeSeries<double> m_means; ///< Time series of means.
    TimeSeries<double> m_moments; ///< Time series of all moments up to MaxMomentOrder.
    std::vector<TimeSeries<ScalarType>> m_results; ///< Interpolated simulation results.
    std::vector<std::string> m_moment_names; ///< Moment names as they are saved in m_moments.
    double m_dt; ///< Interpolation time step.
    std::vector<double> m_sim_time; ///< Simulation time per run.
    double m_advance_time; ///< Time the advance function took in total.
    MomentArray<static_cast<size_t>(Status::Count), regions, MaxMomentOrder>
        m_mom_array{}; // Moment array used to calculate moments and their names.
};

} // namespace smm
} // namespace mio

#endif // MIO_SMM_SIMULATION_SET_H
