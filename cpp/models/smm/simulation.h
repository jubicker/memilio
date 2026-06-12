/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: René Schmieding, Julia Bicker
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

#ifndef MIO_SMM_SIMULATION_H
#define MIO_SMM_SIMULATION_H

#include "memilio/config.h"
#include "memilio/epidemiology/adoption_rate.h"
#include "memilio/utils/time_series.h"
#include "smm/model.h"
#include "smm/parameters.h"
#include "memilio/compartments/simulation.h"
#include <algorithm>
#include <cstddef>
#include <iostream>
#include <map>

namespace mio
{

namespace smm
{

/**
 * @brief A specialized Simulation for mio::smm::Model.
 * @tparam regions The number of regions.
 * @tparam Status An infection state enum.
 */
template <typename FP, size_t regions, class Status>
class Simulation
{

public:
    using Model = smm::Model<FP, regions, Status>;

    /**
     * @brief Set up the simulation for a Stochastic Metapopulation Model.
     * @param[in] model An instance of mio::smm::Model.
     * @param[in] t0 Start time.
     * @param[in] dt Initial Step size.
     */
    Simulation(Model const& model, FP t0 = 0., FP dt = 1.)
        : num_transitions(
              1, Eigen::Matrix<size_t, static_cast<size_t>(Status::Count) * regions,
                               static_cast<size_t>(Status::Count) * regions>::Zero(static_cast<size_t>(Status::Count) *
                                                                                       regions,
                                                                                   static_cast<size_t>(Status::Count) *
                                                                                       regions))
        , num_events(t0, Eigen::VectorXd::Zero(1))
        , m_dt(dt)
        , m_model(std::make_unique<Model>(model))
        , m_result_interpolated(t0, m_model->get_initial_values())
        , m_internal_time(adoption_rates().size() + transition_rates().size(), t0)
        , m_tp_next_event(adoption_rates().size() + transition_rates().size(), t0)
        , m_waiting_times(adoption_rates().size() + transition_rates().size(), 0)
        , m_current_rates(adoption_rates().size() + transition_rates().size(), 0)

    {
        assert(dt > 0);
        assert(m_waiting_times.size() > 0);
        assert(std::all_of(adoption_rates().begin(), adoption_rates().end(), [](auto&& r) {
            return static_cast<size_t>(r.region) < regions;
        }));
        assert(std::all_of(transition_rates().begin(), transition_rates().end(), [](auto&& r) {
            return static_cast<size_t>(r.from) < regions && static_cast<size_t>(r.to) < regions;
        }));
        // initialize (internal) next event times by random values
        for (size_t i = 0; i < m_tp_next_event.size(); i++) {
            m_tp_next_event[i] += mio::ExponentialDistribution<FP>::get_instance()(m_model->get_rng(), 1.0);
        }
    }

    Simulation(const Simulation& other)
        : num_transitions(other.num_transitions)
        , num_events(other.num_events)
        , m_dt(other.m_dt)
        , m_model(std::make_unique<Model>(*other.m_model))
        , m_result_interpolated(other.m_result_interpolated)
        , m_internal_time(other.m_internal_time)
        , m_tp_next_event(other.m_tp_next_event)
        , m_waiting_times(other.m_waiting_times)
        , m_current_rates(other.m_current_rates)
    {
    }

    /**
     * @brief Advance simulation to tmax.
     * This function performs a Gillespie algorithm.
     * @param tmax Next stopping point of simulation.
     */
    Eigen::Ref<Eigen::VectorX<FP>> advance(FP tmax)
    {
        update_current_rates_and_waiting_times();
        size_t next_event = determine_next_event(); // index of the next event
        FP current_time   = m_result_interpolated.get_last_time();
        // set in the past to add a new time point immediately
        FP next_result_time = m_result_interpolated.get_last_time() + m_dt;
        // iterate over time
        while (current_time + m_waiting_times[next_event] < tmax) {
            // update time
            current_time += m_waiting_times[next_event];
            // Regularly save current state
            if (current_time > next_result_time) {
                if (count_transitions) {
                    if (int(next_result_time) > int(num_transitions.size()) - 1) {
                        num_transitions.push_back(
                            Eigen::Matrix<size_t, static_cast<size_t>(Status::Count) * regions,
                                          static_cast<size_t>(Status::Count) *
                                              regions>::Zero(static_cast<size_t>(Status::Count) * regions,
                                                             static_cast<size_t>(Status::Count) * regions));
                    }
                }
                while (current_time > next_result_time) {
                    if (count_events) {
                        num_events.add_time_point(next_result_time, Eigen::VectorXd::Zero(1));
                    }
                    m_result_interpolated.add_time_point(next_result_time);
                    // copy from the previous last value
                    m_result_interpolated.get_last_value() = m_model->populations.get_compartments();
                    next_result_time += m_dt;
                }
            }
            // decide event type by index and perform it
            if (next_event < adoption_rates().size()) {
                // perform adoption event
                const auto& rate = adoption_rates()[next_event];
                m_model->populations[{rate.region, rate.from}] -= 1.0;
                m_model->populations[{rate.region, rate.to}] += 1.0;
                if (count_events) {
                    num_events.get_last_value()[0] += 1;
                }
            }
            else {
                // perform transition event
                const auto& rate = transition_rates()[next_event - adoption_rates().size()];
                m_model->populations[{rate.from, rate.status}] -= 1.0;
                m_model->populations[{rate.to, rate.status}] += 1.0;

                if (count_transitions) {
                    num_transitions.back()(static_cast<size_t>(rate.from) * static_cast<size_t>(Status::Count) +
                                               static_cast<size_t>(rate.status),
                                           static_cast<size_t>(rate.to) * static_cast<size_t>(Status::Count) +
                                               static_cast<size_t>(rate.status)) += 1;
                }
                if (count_events) {
                    num_events.get_last_value()[0] += 1;
                }
            }
            // update internal times
            for (size_t i = 0; i < m_internal_time.size(); i++) {
                m_internal_time[i] += m_current_rates[i] * m_waiting_times[next_event];
            }
            // draw new "next event" time for the occured event
            m_tp_next_event[next_event] += mio::ExponentialDistribution<FP>::get_instance()(m_model->get_rng(), 1.0);
            // precalculate next event
            if (next_event < adoption_rates().size()) {
                update_current_rates_and_waiting_times(adoption_rates()[next_event]);
            }
            else {
                update_current_rates_and_waiting_times(transition_rates()[next_event - adoption_rates().size()]);
            }
            //update_current_rates_and_waiting_times();
            next_event = determine_next_event();
        }
        // copy last result, if no event occurs between last_result_time and tmax
        if (m_result_interpolated.get_last_time() < tmax) {
            while (tmax >=
                   (next_result_time -
                    1e-10)) { // add time points until tmax is reached, with a small tolerance to avoid numerical issues
                m_result_interpolated.add_time_point(std::min(next_result_time, tmax));
                m_result_interpolated.get_last_value() = m_model->populations.get_compartments();
                if (count_transitions) {
                    if (int(next_result_time) > int(num_transitions.size()) - 1) {
                        num_transitions.push_back(
                            Eigen::Matrix<size_t, static_cast<size_t>(Status::Count) * regions,
                                          static_cast<size_t>(Status::Count) *
                                              regions>::Zero(static_cast<size_t>(Status::Count) * regions,
                                                             static_cast<size_t>(Status::Count) * regions));
                    }
                }
                next_result_time += m_dt;
            }
            // update internal times
            for (size_t i = 0; i < m_internal_time.size(); i++) {
                m_internal_time[i] += m_current_rates[i] * (tmax - current_time);
            }
        }
        return m_result_interpolated.get_last_value();
    }

    /**
     * @brief Returns the final simulation result.
     * @return A TimeSeries to represent the final simulation result.
     */
    TimeSeries<FP>& get_result()
    {
        return m_result_interpolated;
    }
    const TimeSeries<FP>& get_result() const
    {
        return m_result_interpolated;
    }

    /**
     * @brief Returns the model used in the simulation.
     */
    const Model& get_model() const
    {
        return *m_model;
    }
    Model& get_model()
    {
        return *m_model;
    }

    std::vector<Eigen::Matrix<size_t, static_cast<size_t>(Status::Count) * regions,
                              static_cast<size_t>(Status::Count) * regions>>
        num_transitions; ///< Matrix counting the number of spatial transitions between regions for each infection state per time step.
    mio::TimeSeries<double>
        num_events; ///< Time series counting the number of events that happened until the next time series time point.

private:
    /**
     * @brief Returns the model's transition rates.
     */
    inline constexpr const typename smm::TransitionRates<FP, Status>::Type& transition_rates()
    {
        return m_model->parameters.template get<smm::TransitionRates<FP, Status>>();
    }

    /**
     * @brief Returns the model's adoption rates.
     */
    inline constexpr const typename smm::AdoptionRates<FP, Status>::Type& adoption_rates()
    {
        return m_model->parameters.template get<smm::AdoptionRates<FP, Status>>();
    }

    /**
     * @brief Calculate current values for m_current_rates and m_waiting_times.
     */
    inline void update_current_rates_and_waiting_times()
    {
        size_t i                  = 0; // shared index for iterating both rates
        double seasonality_factor = 1.0;
        if (m_model->parameters.template get<SeasonalityRho>().size() >= 1) {
            seasonality_factor = calculate_seasonality_factor(m_result_interpolated.get_last_time());
        }
        for (const auto& rate : adoption_rates()) {
            m_current_rates[i] = m_model->evaluate(rate, m_model->populations.get_compartments(), seasonality_factor);
            m_waiting_times[i] = (m_current_rates[i] > 0)
                                     ? (m_tp_next_event[i] - m_internal_time[i]) / m_current_rates[i]
                                     : std::numeric_limits<FP>::max();
            i++;
        }
        for (const auto& rate : transition_rates()) {
            m_current_rates[i] = m_model->evaluate(rate, m_model->populations.get_compartments());
            m_waiting_times[i] = (m_current_rates[i] > 0)
                                     ? (m_tp_next_event[i] - m_internal_time[i]) / m_current_rates[i]
                                     : std::numeric_limits<FP>::max();
            i++;
        }
    }

    inline void update_adoption_rate(double seasonality_factor, const AdoptionRate<FP, Status>& rate, size_t rate_index)
    {
        m_current_rates[rate_index] =
            m_model->evaluate(rate, m_model->populations.get_compartments(), seasonality_factor);
    }

    inline void update_transition_rate(const TransitionRate<FP, Status>& rate, size_t rate_index)
    {
        m_current_rates[rate_index] = m_model->evaluate(rate, m_model->populations.get_compartments());
    }

    /**
     * @brief Update current values for m_current_rates and m_waiting_times given the last adoption event.
     */
    inline void update_current_rates_and_waiting_times(const AdoptionRate<FP, Status>& last_event_rate)
    {
        double seasonality_factor = 1.0;
        if (m_model->parameters.template get<SeasonalityRho>().size() >= 1) {
            seasonality_factor = calculate_seasonality_factor(m_result_interpolated.get_last_time());
        }
        size_t i = 0; // shared index for iterating both rates
        for (const auto& rate : adoption_rates()) {
            if (rate.influences.size() == 0) { // First-order adoptions
                if (rate.region ==
                    last_event_rate.region) // First-order adoption rates only change populations in within one region
                {
                    if (rate.from == last_event_rate.from || rate.from == last_event_rate.to) {
                        update_adoption_rate(seasonality_factor, rate, i);
                    }
                }
            }
            else { // Second-order adoption rates also have to be recalculated if one of their influences has changed
                if ((rate.region == last_event_rate.region &&
                     (rate.from == last_event_rate.from || rate.from == last_event_rate.to)) ||
                    (std::find_if(rate.influences.begin(), rate.influences.end(),
                                  [last_event_rate](const Influence<FP, Status>& influence) {
                                      return (influence.status == last_event_rate.from ||
                                              influence.status == last_event_rate.to) &&
                                             (influence.region == last_event_rate.region);
                                  }) != rate.influences.end())) {
                    update_adoption_rate(seasonality_factor, rate, i);
                }
            }
            m_waiting_times[i] = (m_current_rates[i] > 0)
                                     ? (m_tp_next_event[i] - m_internal_time[i]) / m_current_rates[i]
                                     : std::numeric_limits<FP>::max();
            i++;
        }
        for (const auto& rate : transition_rates()) {
            if ((rate.status == last_event_rate.from || rate.status == last_event_rate.to) &&
                (rate.from == last_event_rate.region)) {
                update_transition_rate(rate, i);
            }
            m_waiting_times[i] = (m_current_rates[i] > 0)
                                     ? (m_tp_next_event[i] - m_internal_time[i]) / m_current_rates[i]
                                     : std::numeric_limits<FP>::max();
            i++;
        }
    }

    /**
     * @brief Update current values for m_current_rates and m_waiting_times given the last transition event.
     */
    inline void update_current_rates_and_waiting_times(const TransitionRate<FP, Status>& last_event_rate)
    {
        double seasonality_factor = 1.0;
        if (m_model->parameters.template get<SeasonalityRho>().size() >= 1) {
            seasonality_factor = calculate_seasonality_factor(m_result_interpolated.get_last_time());
        }
        size_t i = 0; // shared index for iterating both rates
        for (const auto& rate : adoption_rates()) {
            if (rate.influences.size() == 0) { // First-order adoptions
                if (rate.region == last_event_rate.from ||
                    rate.region == last_event_rate.to) // Check if first-order adoption lies within midoefies region
                {
                    if (rate.from == last_event_rate.status) {
                        update_adoption_rate(seasonality_factor, rate, i);
                    }
                }
            }
            else { // Second-order adoption rates also have to be recalculated if one of their influences has changed
                if ((rate.from == last_event_rate.status &&
                     (rate.region == last_event_rate.from || rate.region == last_event_rate.to)) ||
                    (std::find_if(rate.influences.begin(), rate.influences.end(),
                                  [last_event_rate](const Influence<FP, Status>& influence) {
                                      return (influence.region == last_event_rate.from ||
                                              influence.region == last_event_rate.to) &&
                                             (influence.status == last_event_rate.status);
                                  }) != rate.influences.end())) {
                    update_adoption_rate(seasonality_factor, rate, i);
                }
            }
            m_waiting_times[i] = (m_current_rates[i] > 0)
                                     ? (m_tp_next_event[i] - m_internal_time[i]) / m_current_rates[i]
                                     : std::numeric_limits<FP>::max();
            i++;
        }
        for (const auto& rate : transition_rates()) {
            if ((rate.status == last_event_rate.status) &&
                (rate.from == last_event_rate.from || rate.from == last_event_rate.to)) {
                update_transition_rate(rate, i);
            }
            m_waiting_times[i] = (m_current_rates[i] > 0)
                                     ? (m_tp_next_event[i] - m_internal_time[i]) / m_current_rates[i]
                                     : std::numeric_limits<FP>::max();
            i++;
        }
    }

    /**
     * @brief Get next event i.e. event with the smallest waiting time.
     */
    inline size_t determine_next_event()
    {
        return std::distance(m_waiting_times.begin(), std::min_element(m_waiting_times.begin(), m_waiting_times.end()));
    }

    double gaussian(double t, double t_peak_season, int season)
    {
        double sigma = m_model->parameters.template get<SeasonalitySigma>()[season];
        return 1. / (sigma * sigma * std::sqrt(2 * M_PI)) * std::exp(-0.5 * std::pow((t - t_peak_season) / sigma, 2));
    }

    double zeta(double t, int season)
    {
        double tau_minus = m_model->parameters.template get<FirstSeasonStartDay>() -
                           m_model->parameters.template get<StartDay>() + season * 365.;
        double tau_plus = tau_minus + 365.;
        if (t <= tau_minus + 30) {
            return (t - tau_minus + 30.) / 60.;
        }
        else if (t >= tau_plus - 30. && t <= tau_plus) {
            return (-t + tau_plus + 30.) / 60.;
        }
        return 1;
    }

    double calculate_seasonality_factor(double t)
    {
        int season           = int((t + m_model->parameters.template get<StartDay>() -
                          m_model->parameters.template get<FirstSeasonStartDay>()) /
                         365.);
        double t_peak_season = season * 365. + m_model->parameters.template get<SeasonalityPeak>()[season] -
                               m_model->parameters.template get<StartDay>();
        double rho = m_model->parameters.template get<SeasonalityRho>()[season];
        return rho + (1 - rho) * gaussian(t, t_peak_season, season) / gaussian(0, 0, season) * zeta(t, season);
    }

    FP m_dt; ///< Initial step size
    std::unique_ptr<Model> m_model; ///< Pointer to the model used in the simulation.
    mio::TimeSeries<FP> m_result_interpolated; ///< Interpolated result.
    std::vector<FP> m_internal_time; ///< Internal times of all poisson processes (aka T_k).
    std::vector<FP> m_tp_next_event; ///< Internal time points of next event i after m_internal[i] (aka P_k).
    std::vector<FP> m_waiting_times; ///< External times between m_internal_time and m_tp_next_event.
    std::vector<FP> m_current_rates; ///< Current values of both types of rates i.e. adoption and transition rates.

    const bool count_transitions = true;
    const bool count_events      = true;
};

} //namespace smm
} // namespace mio

#endif
