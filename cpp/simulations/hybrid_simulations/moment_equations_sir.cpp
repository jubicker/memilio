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

#include "memilio/compartments/simulation.h"
#include "memilio/config.h"
#include "memilio/data/analyze_result.h"
#include "memilio/math/integrator.h"
#include "memilio/utils/time_series.h"
#include "moment_equation_array.h"
#include "models/ode_sir/infection_state.h"
#include "models/ode_sir/model.h"
#include "setup.cpp"
#include <boost/math/special_functions/math_fwd.hpp>
#include <cassert>
#include <vector>
#include <cstddef>

#include <boost/math/special_functions/binomial.hpp>

template double boost::math::binomial_coefficient<double>(unsigned, unsigned);

namespace SIR
{
template <int Order>
class Moments
{
public:
    Moments(std::vector<double> initial_populations)
        : values(initial_populations)
    {
    }

    void get_derivatives(Eigen::Ref<const Eigen::VectorX<ScalarType>> y, ScalarType /*t*/,
                         Eigen::Ref<Eigen::VectorX<ScalarType>> dydt, int closure_order) const
    {
        size_t mu_S = static_cast<size_t>(mio::osir::InfectionState::Susceptible);
        size_t mu_I = static_cast<size_t>(mio::osir::InfectionState::Infected);
        size_t mu_R = static_cast<size_t>(mio::osir::InfectionState::Recovered);
        assert(closure_order >= 2);
        double M_110 = closure_order > 2 ? y[values.get_moment_index({MyIndex(1), MyIndex(1), MyIndex(0)})] : 0;
        dydt[mu_S]   = -setup::lambda * y[mu_S] * y[mu_I] - setup::lambda * M_110;
        dydt[mu_I]   = setup::lambda * y[mu_S] * y[mu_I] - setup::gamma * y[mu_I] + setup::lambda * M_110;
        dydt[mu_R]   = setup::gamma * y[mu_I];
        if (closure_order > 2) {
            for (int i_S = 0; i_S <= closure_order; i_S++) {
                for (int i_I = 0; i_I <= closure_order; i_I++) {
                    for (int i_R = 0; i_R <= closure_order; i_R++) {
                        if (i_S + i_I + i_R == 0 || i_S + i_I + i_R == 1) {
                            continue;
                        }
                        else {
                            size_t index = values.get_moment_index({MyIndex(i_S), MyIndex(i_I), MyIndex(i_R)});
                            dydt[index]  = 0;
                            for (int l_S = 0; l_S <= i_S; l_S++) {
                                for (int l_I = 0; l_I <= i_I; l_I++) {
                                    if (!(l_S + l_I + i_R == 1)) {
                                        double M_lS_lI_iR =
                                            y[values.get_moment_index({MyIndex(l_S), MyIndex(l_I), MyIndex(i_R)})];
                                        if (l_S + l_I + i_R == 0) {
                                            M_lS_lI_iR = 1.;
                                        }
                                        dydt[index] += setup::lambda * y[mu_S] * y[mu_I] *
                                                       boost::math::binomial_coefficient<double>(i_S, l_S) *
                                                       boost::math::binomial_coefficient<double>(i_I, l_I) *
                                                       std::pow(-1, i_S - l_S) * M_lS_lI_iR;
                                    }
                                    if (!(l_S + 1 + l_I + i_R == 1)) {
                                        double M_lSp1_lI_iR =
                                            y[values.get_moment_index({MyIndex(l_S + 1), MyIndex(l_I), MyIndex(i_R)})];
                                        if (l_S + 1 + l_I + i_R == 0) {
                                            M_lSp1_lI_iR = 1.;
                                        }
                                        dydt[index] += setup::lambda * y[mu_I] *
                                                       boost::math::binomial_coefficient<double>(i_S, l_S) *
                                                       boost::math::binomial_coefficient<double>(i_I, l_I) *
                                                       std::pow(-1, i_S - l_S) * M_lSp1_lI_iR;
                                    }
                                    if (!(l_S + l_I + 1 + i_R == 1)) {
                                        double M_lS_lIp1_iR =
                                            y[values.get_moment_index({MyIndex(l_S), MyIndex(l_I + 1), MyIndex(i_R)})];
                                        if (l_S + l_I + 1 + i_R == 0) {
                                            M_lS_lIp1_iR = 1.;
                                        }
                                        dydt[index] += setup::lambda * y[mu_S] *
                                                       boost::math::binomial_coefficient<double>(i_S, l_S) *
                                                       boost::math::binomial_coefficient<double>(i_I, l_I) *
                                                       std::pow(-1, i_S - l_S) * M_lS_lIp1_iR;
                                    }
                                    if (!(l_S + 1 + l_I + 1 + i_R == 1)) {
                                        double M_lSp1_lIp1_iR = y[values.get_moment_index(
                                            {MyIndex(l_S + 1), MyIndex(l_I + 1), MyIndex(i_R)})];
                                        if (l_S + 1 + l_I + 1 + i_R == 0) {
                                            M_lSp1_lIp1_iR = 1.;
                                        }
                                        dydt[index] += setup::lambda *
                                                       boost::math::binomial_coefficient<double>(i_S, l_S) *
                                                       boost::math::binomial_coefficient<double>(i_I, l_I) *
                                                       std::pow(-1, i_S - l_S) * M_lSp1_lIp1_iR;
                                    }
                                }
                            }
                            for (int l_I = 0; l_I <= i_I; l_I++) {
                                for (int l_R = 0; l_R <= i_R; l_R++) {
                                    if (!(i_S + l_I + l_R == 1)) {
                                        double M_iS_lI_lR =
                                            y[values.get_moment_index({MyIndex(i_S), MyIndex(l_I), MyIndex(l_R)})];
                                        if (i_S + l_I + l_R == 0) {
                                            M_iS_lI_lR = 1.;
                                        }
                                        dydt[index] += setup::gamma * y[mu_I] *
                                                       boost::math::binomial_coefficient<double>(i_I, l_I) *
                                                       boost::math::binomial_coefficient<double>(i_R, l_R) *
                                                       std::pow(-1, i_I - l_I) * M_iS_lI_lR;
                                    }

                                    if (!(i_S + l_I + 1 + l_R == 1)) {
                                        double M_iS_lIp1_lR =
                                            y[values.get_moment_index({MyIndex(i_S), MyIndex(l_I + 1), MyIndex(l_R)})];
                                        if (i_S + l_I + 1 + l_R == 0) {
                                            M_iS_lIp1_lR = 1.;
                                        }
                                        dydt[index] += setup::gamma *
                                                       boost::math::binomial_coefficient<double>(i_I, l_I) *
                                                       boost::math::binomial_coefficient<double>(i_R, l_R) *
                                                       std::pow(-1, i_I - l_I) * M_iS_lIp1_lR;
                                    }
                                }
                            }
                            if (!(i_S - 1 + i_I + i_R == 1 || i_S - 1 + i_I + i_R < 0)) {
                                double M_iSm1_iI_iR =
                                    y[values.get_moment_index({MyIndex(i_S - 1), MyIndex(i_I), MyIndex(i_R)})];
                                double M_iS_iIm1_iR =
                                    y[values.get_moment_index({MyIndex(i_S), MyIndex(i_I - 1), MyIndex(i_R)})];
                                double M_iS_iI_iRm1 =
                                    y[values.get_moment_index({MyIndex(i_S), MyIndex(i_I), MyIndex(i_R - 1)})];
                                if (i_S - 1 + i_I + i_R == 0) {
                                    M_iSm1_iI_iR = 1.;
                                    M_iS_iIm1_iR = 1.;
                                    M_iS_iI_iRm1 = 1.;
                                }
                                dydt[index] -= i_S * dydt[mu_S] * M_iSm1_iI_iR;
                                dydt[index] -= i_I * dydt[mu_I] * M_iS_iIm1_iR;
                                dydt[index] -= i_R * dydt[mu_R] * M_iS_iI_iRm1;
                            }
                        }
                    }
                }
            }
        }
    }

    MomentEquationArray<mio::osir::InfectionState, Order> values;
};
} // namespace SIR

int main()
{
    const size_t order = 2;
    SIR::Moments<order> moments({setup::total_population - setup::I0, setup::I0, 0});
    mio::TimeSeries<ScalarType> result(setup::t0, moments.values.get_values());
    mio::OdeIntegrator<ScalarType> integrator(std::make_shared<mio::DefaultIntegratorCore<ScalarType>>());

    integrator.advance(
        [order, moments](auto&& y, auto&& t, auto&& dydt) {
            moments.get_derivatives(y, t, dydt, order);
        },
        setup::tmax, setup::dt, result);

    std::string output_file = "/home/bick_ju/Documents/TemporalHybridModel/sim_outputs/moments/sir_order2.csv";

    int num_steps = static_cast<int>(setup::tmax / setup::dt) + 1;

    std::vector<double> interpolation_tps(num_steps);

    for (int i = 0; i < num_steps; ++i) {
        interpolation_tps[i] = i * 0.1;
    }
    const std::vector<std::string> second_order = {"muS", "muI", "muR", "M_000", "M_100", "M_010", "M_001"};
    auto finished = mio::interpolate_simulation_result(result, interpolation_tps).export_csv(output_file, second_order);
    return 0;
}
