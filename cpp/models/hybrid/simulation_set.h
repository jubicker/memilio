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

#ifndef MIO_HYBRID_SIMULATION_SET_H
#define MIO_HYBRID_SIMULATION_SET_H

#include "hybrid/temporal_hybrid_model.h"
#include "memilio/config.h"
#include "memilio/data/analyze_result.h"
#include "memilio/timer/basic_timer.h"
#include "memilio/utils/time_series.h"
#include "smm/simulation.h"
#include "smm_moments/simulation.h"
#include <cstddef>
#include <vector>
namespace mio
{
namespace hybrid
{

/**
 * @brief This class performs multiple simulations of the hybrid SMM-ODE model, see hybrid::TemporalHybridSimulation.
 * @tparam regions The number of regions.
 * @tparam Status The used infection state enum.
 * @tparam MaxMomentOrder The maximum order up to which moments are calculated as result.
 */
template <size_t regions, class Status, size_t MaxMomentOrder>
class SimulationSet
{
public:
    using Simulation =
        hybrid::TemporalHybridSimulation<smm::Simulation<ScalarType, regions, Status>,
                                         smm_moments::Simulation<regions, 2>, TimeSeries<double>, TimeSeries<double>>;
    using Model1 = smm::Model<ScalarType, regions, Status>;
    using Model2 = smm_moments::Model<regions, 2>;
    using result1_function =
        std::function<TimeSeries<double>(smm::Simulation<ScalarType, regions, Status>&, ScalarType t)>;
    using result2_function    = std::function<TimeSeries<double>(smm_moments::Simulation<regions, 2>&, ScalarType t)>;
    using switching_condition = Simulation::switching_condition;

    /**
     * @brief Constructor for SimulationSet.
     * @param num_runs The number of simulations to run.
     * @param sim The simulation to run.
     * @param dt The time step for interpolation.
     */
    SimulationSet(size_t num_runs, const Model1& model1, const Model2& model2, const result1_function& result1,
                  const result2_function& result2, ScalarType t0, ScalarType dt, ScalarType dt_switch)
        : m_model1s(num_runs, model1)
        , m_model2s(num_runs, model2)
        , m_results(num_runs, TimeSeries<ScalarType>(static_cast<size_t>(Status::Count) * regions))
        , m_means(static_cast<size_t>(Status::Count) * regions)
        , m_moments(1)
        , m_sim_time(num_runs, 0.0)
        , m_dt(dt)
        , m_t(t0)
    {
        m_moments       = TimeSeries<double>(m_mom_array.names_up_to_order(MaxMomentOrder).size());
        m_moments_names = m_mom_array.names_up_to_order(MaxMomentOrder);

        u_int32_t seed = 0;
        for (size_t run = 0; run < num_runs; ++run) {
            m_model1s[run].get_rng().seed({seed});
            auto sim1 = smm::Simulation<ScalarType, regions, Status>(m_model1s[run], t0, dt);
            auto sim2 = smm_moments::Simulation<regions, 2>(m_model2s[run], t0, dt);
            m_simulations.push_back(
                Simulation(std::move(sim1), std::move(sim2), result1, result2, true, t0, dt_switch));
            seed++;
        }
    }

    void advance(double tmax, const switching_condition& switch_model)
    {
        timing::BasicTimer total_timer;
        total_timer.start();
// Run simulations
#pragma omp parallel for
        for (size_t run = 0; run < m_simulations.size(); ++run) {
            timing::BasicTimer timer;
            timer.start();
            m_simulations[run].advance(tmax, switch_model);
            timer.stop();
            m_sim_time[run] += timer.get_elapsed_time();
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

        // Get last output time
        if (m_means.get_num_time_points() > 0) {
            t = m_means.get_last_time();
        }

        // Interpolation time steps
        int num_steps = static_cast<int>((m_t - t) / m_dt);
        if (num_steps > 0) {
            num_steps += 1;

            std::vector<double> interpolation_tps(num_steps);
            for (int i = 0; i < num_steps; ++i) {
                interpolation_tps[i] = t + i * m_dt;
            }

            for (size_t run = 0; run < m_simulations.size(); ++run) {
                // Interpolate SMM and ODE results separately
                auto interpolated_smm_result =
                    interpolate_smm_simulation_result(m_simulations[run].get_result_model1(), interpolation_tps);
                auto interpolated_ode_result =
                    interpolate_simulation_result(m_simulations[run].get_result_model2(), interpolation_tps);
                auto sim_ts     = m_simulations[run].get_result_model2();
                auto num_points = static_cast<size_t>(sim_ts.get_num_time_points());
                for (size_t i = 0; i < num_points; i++) {
                    printf("\n%.14f ", sim_ts.get_time(i));
                    Eigen::VectorX<ScalarType> res_j = sim_ts.get_value(i);
                    for (size_t j = 0; j < (size_t)res_j.size(); j++) {
                        printf(" %.14f", res_j[j]);
                    }
                }
                auto merged_ts = merge_time_series(interpolated_smm_result, interpolated_ode_result).value();
                while (merged_ts.get_num_time_points() > 1) {
                    if (m_results[run].get_num_time_points() == 0 ||
                        m_results[run].get_last_time() < merged_ts.get_time(0)) {
                        m_results[run].add_time_point(merged_ts.get_time(0), merged_ts.get_value(0));
                    }
                    merged_ts.remove_time_point(0);
                }
                m_results[run].add_time_point(merged_ts.get_last_time(), merged_ts.get_last_value());
            }

            //Update means and moments
            update_means();
            update_moments();
        }
    }

    /**
     * @brief Get all simulations.
     */
    std::vector<Simulation>& get_simulations()
    {
        return m_simulations;
    }
    const std::vector<Simulation>& get_simulations() const
    {
        return m_simulations;
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
        return m_moments_names;
    }
    const std::vector<std::string> get_moment_names() const
    {
        return m_moments_names;
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

private:
    /**
     * @brief Calculate mean time series from simulation results.
     */
    void update_means()
    {
        const size_t num_elements = static_cast<size_t>(Status::Count) * regions;
        for (auto t = m_t_index; t < m_results[0].get_num_time_points(); ++t) {
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
        Eigen::Matrix<ScalarType, 1, static_cast<size_t>(Status::Count) * regions> means = values.colwise().mean();
        double moment                                                                    = 0.0;
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
    void update_moments()
    {
        const size_t num_elements = static_cast<size_t>(Status::Count) * regions;
        for (int t = m_t_index; t < m_results[0].get_num_time_points(); ++t) {
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
            m_t_index += 1;
        }
    }

    std::vector<Model1> m_model1s; ///< SMM models for all runs.
    std::vector<Model2> m_model2s; ///< Moment models for all runs.
    std::vector<Simulation> m_simulations; ///< Hybrid simulations.
    std::vector<TimeSeries<double>> m_results; ///< Interpolated simulation results.
    TimeSeries<double> m_means; ///< Time series of means.
    TimeSeries<double> m_moments; ///< Time series of all moments up to MaxMomentOrder.
    std::vector<std::string> m_moments_names; ///< Moment names as they are saved in m_moments.
    std::vector<double> m_sim_time; ///< Simulation time per run.
    double m_dt; ///< Interpolation time step.
    double m_t; ///< Current time.
    int m_t_index; ///< Current result time point index.
    MomentArray<static_cast<size_t>(Status::Count), regions, MaxMomentOrder>
        m_mom_array{}; // Moment array used to calculate moments and their names.
};

} //namespace hybrid
} //namespace mio

#endif // MIO_HYBRID_SIMULATION_SET_H
