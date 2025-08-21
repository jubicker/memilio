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
#include "memilio/math/euler.h"
#include "memilio/math/integrator.h"
#include "memilio/utils/time_series.h"
#include "moment_equation_array.h"
#include "models/ode_sir/infection_state.h"
#include "models/ode_sir/model.h"
#include "setup.cpp"
#include <boost/math/special_functions/math_fwd.hpp>
#include <cassert>
#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include <cstddef>

int main()
{
    std::string save_file = "/home/bick_ju/Documents/TemporalHybridModel/sim_outputs/moments/I0=" +
                            std::to_string(static_cast<int>(setup::I0)) + "RKnormal/";
    setup::save_setup(save_file + "setup.csv");

    const size_t closing_order = 3;
    SIR::Moments<closing_order> moments({setup::total_population - setup::I0, setup::I0, 0});
    mio::TimeSeries<ScalarType> result(setup::t0, moments.values.get_values());

    //RK54 integrator
    mio::OdeIntegrator<ScalarType> integrator(std::make_shared<mio::DefaultIntegratorCore<ScalarType>>(
        1e-10, 1e-5, std::numeric_limits<double>::min(),
        std::numeric_limits<double>::max())); //abs_tol, rel_tol, dt_min, dt_max //setup::dt / 1000.

    //mio::OdeIntegrator<ScalarType> integrator(std::make_shared<mio::EulerIntegratorCore<ScalarType>>());

    double dt_integrator = setup::dt;
    integrator.advance(
        [closing_order, moments](auto&& y, auto&& t, auto&& dydt) {
            moments.get_derivatives(y, t, dydt, closing_order, setup::lambda, setup::gamma);
        },
        setup::tmax, dt_integrator, result);

    std::string output_file = save_file + "sir_order" + std::to_string(closing_order) + ".csv";

    auto result_moments = moments.get_moment_ts_and_names(result, closing_order);
    auto finished       = result_moments.first.export_csv(output_file, result_moments.second);

    return 0;
}
