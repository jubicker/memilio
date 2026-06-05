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
#include "memilio/timer/auto_timer.h"
#include "memilio/timer/basic_timer.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/random_number_generator.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/mioomp.h"
#include "ode_sir/infection_state.h"
#include "smm/simulation.h"
#include "simulations/hybrid_simulations/sir_metapop/library/moment_array.h"
#include "memilio/data/analyze_result.h"
#include <Eigen/src/Core/Matrix.h>
#include <cstddef>
#include <limits>
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
    using MatrixType = Eigen::Matrix<ScalarType, Eigen::Dynamic, static_cast<size_t>(Status::Count) * regions>;

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
        calculate_outputs();
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
            m_sims[run].advance(tmax);
        }

        m_t = tmax;

        calculate_outputs();
    }

    /**
     * @brief Calculate interpolated results, mean and moment time series. Only new time points since last call are processed, so this can be called multiple times during a simulation to get intermediate results without reprocessing old time points.
     */
    void calculate_outputs()
    {
        // Struct has time and compartment values of one simulation for the time point
        struct TPValues {
            double time;
            Eigen::VectorXd value;
        };
        // For each run: Vector with new TPValues since last calculate_outputs() call (i.e. since last interpolation)
        std::vector<std::vector<TPValues>> new_points(m_sims.size());
#ifdef MEMILIO_ENABLE_OPENMP
#pragma omp parallel for
#endif
        for (size_t run = 0; run < m_sims.size(); ++run) {
            auto& sim_ts          = m_sims[run].get_result();
            const size_t n_points = sim_ts.get_num_time_points();
            if (n_points == 0) {
                continue;
            }
            new_points[run].reserve(n_points);
            // Copy simulation time points to small vector
            for (size_t i = 0; i < n_points; ++i) {
                double time_i = sim_ts.get_time(i);
                auto v        = sim_ts.get_value(i);
                if (m_results[run].get_num_time_points() != 0 &&
                    m_results[run].get_last_time() >=
                        time_i) { // Only add new time points since last calculate_outputs() call, so we don't have to recalculate means and moments for old time points
                    continue;
                }
                new_points[run].push_back({time_i, Eigen::VectorXd::Map(v.data(), v.size())});
            }
            // append all new points to m_results[run]
            for (const auto& tp_value : new_points[run]) {
                if (m_results[run].get_num_time_points() == 0 || m_results[run].get_last_time() < tp_value.time) {
                    m_results[run].add_time_point(tp_value.time, tp_value.value);
                }
            }
            // keep only last point in sim_ts to avoid reprocessing
            const auto& lastp = new_points[run].back();
            sim_ts            = mio::TimeSeries<double>(sim_ts.get_num_elements());
            sim_ts.add_time_point(lastp.time, lastp.value);
        }

        // Prepare sizes
        const int n_runs = static_cast<int>(m_sims.size());
        const int n_cols = last_results.cols();

        // Get moment indices and initialize vector with moments that do not have to be recalculated
        auto& moment_indices = m_mom_array.get_indices();
        std::vector<char> skip_moment(moment_indices.size(), 0);

        // Iterate over all new time points
        for (size_t tp_index = 0; tp_index < new_points[0].size(); ++tp_index) {
            double time_k = new_points[0][tp_index].time;
#ifdef MEMILIO_ENABLE_OPENMP
#pragma omp parallel
#endif
            { // One parallel region per time step
#ifdef MEMILIO_ENABLE_OPENMP
                const int tid      = omp_get_thread_num();
                const int nthreads = omp_get_num_threads();
#else
                const int tid      = 0;
                const int nthreads = 1;
#endif

                // Each thread updates distinct rows (runs) of last_results.
                for (int run = tid; run < n_runs; run += nthreads) {
                    // copy values into last_results row
                    for (int c = 0; c < n_cols; ++c) {
                        last_results(run, c) = new_points[run][tp_index].value[c];
                    }
                }
                // compute means and centered_values once
#ifdef MEMILIO_ENABLE_OPENMP
#pragma omp barrier
#pragma omp single
#endif
                {
                    means = last_results.colwise().mean();
                    if (m_means.get_num_time_points() == 0 || m_means.get_last_time() < time_k) {
                        m_means.add_time_point(time_k, means);
                    }
                    else {
                        m_means.get_last_value() = means;
                    }
                    // compute centered values
                    centered_values = last_results.rowwise() - means;

                    // detect active columns (any non-zero across runs) to skip empty regions
                    Eigen::RowVectorXd maxabs = last_results.cwiseAbs().colwise().maxCoeff();
                    std::vector<char> active_col(static_cast<size_t>(maxabs.size()));
                    for (int c = 0; c < static_cast<int>(maxabs.size()); ++c) {
                        active_col[static_cast<size_t>(c)] =
                            (maxabs[c] >
                             0.0); // if column is all zero, we can skip it in moment calculations and mark moments that only reference this column as zero
                    }
                    // Rescale cached powers for centered values only if shape/order changed
                    const size_t max_order = MaxMomentOrder;
                    if (m_col_pows.size() != (max_order + 1) || m_col_pows[0].rows() != centered_values.rows() ||
                        m_col_pows[0].cols() != centered_values.cols()) {
                        m_col_pows.assign(static_cast<size_t>(max_order) + 1,
                                          MatrixType::Zero(centered_values.rows(), centered_values.cols()));
                    }
                    // compute integer powers into preallocated matrices (no alloc)
                    m_col_pows[0].setOnes(); // Power 0 is always 1
                    if (max_order >= 1) {
                        m_col_pows[1] = centered_values;
                    }
                    for (size_t p = 2; p <= max_order; ++p) {
                        m_col_pows[p].noalias() = m_col_pows[p - 1].cwiseProduct(centered_values);
                    }
                    // set inactive columns in all powers to zero
                    const int cols = static_cast<int>(m_col_pows[0].cols());
                    for (size_t p = 0; p < m_col_pows.size(); ++p) {
                        for (int c = 0; c < cols; ++c) {
                            if (!active_col[static_cast<size_t>(c)]) {
                                m_col_pows[p].col(c).setZero();
                            }
                        }
                    }

                    // Determine which moments are trivially zero because they reference only inactive columns
                    for (size_t mi = 0; mi < moment_indices.size(); ++mi) {
                        const auto& idx_arr = moment_indices[mi];
                        bool active         = true;
                        for (int c = 0; c < static_cast<int>(idx_arr.size()); ++c) {
                            if (idx_arr[static_cast<size_t>(c)] > 0 && !active_col[static_cast<size_t>(c)]) {
                                active = false;
                                break;
                            }
                        }
                        skip_moment[mi] = !active;
                    }
                }
#ifdef MEMILIO_ENABLE_OPENMP
#pragma omp barrier
#endif

                // Calculate moments in parallel across indices using precomputed powers + skip trivial zeros
#ifdef MEMILIO_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif
                for (int i = 0; i < static_cast<int>(moment_indices.size()); ++i) {
                    if (skip_moment[static_cast<size_t>(i)]) {
                        m_mom_array[moment_indices[static_cast<size_t>(i)]] = 0.0;
                        continue;
                    }
                    auto& idx_arr        = moment_indices[static_cast<size_t>(i)];
                    m_mom_array[idx_arr] = calculate_moment_from_pows(m_col_pows, idx_arr);
                }

                // single thread writes m_moments
#ifdef MEMILIO_ENABLE_OPENMP
#pragma omp barrier
#pragma omp single
#endif
                {
                    auto moment_values   = m_mom_array.get_values();
                    Eigen::VectorXd data = Eigen::VectorXd::Map(moment_values.data(), moment_values.size());
                    if (m_moments.get_num_time_points() == 0 || m_moments.get_last_time() < time_k) {
                        m_moments.add_time_point(time_k, data);
                    }
                    else {
                        m_moments.get_last_value() = data;
                    }
                }
            }
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
        auto& moment_indices = m_mom_array.get_indices();
#ifdef MEMILIO_ENABLE_OPENMP
#pragma omp parallel for
#endif
        for (size_t i = 0; i < moment_indices.size(); ++i) {
            bool is_var = std::count(moment_indices[i].begin(), moment_indices[i].end(), 2) == 1 &&
                          std::count(moment_indices[i].begin(), moment_indices[i].end(), 0) ==
                              static_cast<size_t>(osir::InfectionState::Count) * regions - 1;
            if (!is_var) {
                continue;
            }
            size_t index = std::distance(moment_indices[i].begin(),
                                         std::find(moment_indices[i].begin(), moment_indices[i].end(), 2));
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
        auto& moment_indices = m_mom_array.get_indices();
#ifdef MEMILIO_ENABLE_OPENMP
#pragma omp parallel for
#endif
        for (size_t i = 0; i < moment_indices.size(); ++i) {
            bool is_var = std::count(moment_indices[i].begin(), moment_indices[i].end(), 2) == 1 &&
                          std::count(moment_indices[i].begin(), moment_indices[i].end(), 0) ==
                              static_cast<size_t>(osir::InfectionState::Count) * regions - 1;
            if (!is_var) {
                continue;
            }
            size_t index = std::distance(moment_indices[i].begin(),
                                         std::find(moment_indices[i].begin(), moment_indices[i].end(), 2));
            auto last_tp = m_moments.get_last_time();
            if (m_moments.get_num_time_points() > 1) {
                auto second_last_tp_index = m_moments.get_num_time_points() - 2;
                auto second_last_tp       = m_moments.get_time(second_last_tp_index);
                auto y_second_last        = m_moments.get_value(second_last_tp_index);
                auto y_last               = m_moments.get_last_value();
                vars_gradient[index]      = (y_last[i] - y_second_last[i]) / (last_tp - second_last_tp);
            }
            else {
                vars_gradient[index] = std::numeric_limits<double>::
                    max(); // If there is only one time point, the gradient is set to max double value to ensure that it is above any reasonable threshold for switching conditions
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

    /**
     * @brief Calculates moment defined by indices for the given realizations of a stochastic variable using precomputed powers to avoid repeated pow calls and allocations.
     * @param[in] pows The precomputed column-wise integer powers of the centered values (mean - realization)
     * @param[in] indices The indices defining the moment to be calculated.
     */
    double
    calculate_moment_from_pows(const std::vector<MatrixType>& pows,
                               const std::array<int, static_cast<size_t>(Status::Count) * regions>& indices) const
    {
        const int rows = static_cast<int>(pows[0].rows());
        const int cols = static_cast<int>(pows[0].cols());
        double moment  = 0.0;
        for (int i = 0; i < rows; ++i) {
            double prod = 1.0;
            for (int c = 0; c < cols; ++c) {
                const int p = indices[static_cast<size_t>(c)];
                if (p != 0) {
                    prod *= pows[static_cast<size_t>(p)](i, c);
                }
            }
            moment += prod;
        }
        return moment / static_cast<double>(rows);
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
    std::vector<MatrixType>
        m_col_pows; ///< Cached column powers of centered values used for moment calculation to avoid repeated pow calls and allocations.
};

} // namespace smm
} // namespace mio

#endif // MIO_SMM_SIMULATION_SET_H
