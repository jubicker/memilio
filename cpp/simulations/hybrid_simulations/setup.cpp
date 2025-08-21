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

#include "setup.h"
#include "memilio/utils/logging.h"
#include <string>

namespace setup
{
// Lambda is the rate from S->I (SIR model) or S->E (SEIR/SECIR model)
double lambda = 0.0001;
//gamma = 1 / T_I is the rate from I->R (SIR model)
double gamma = 1. / 5.;

void save_setup(std::string filename)
{
    auto file = fopen(filename.c_str(), "w");
    if (file == NULL) {
        mio::log(mio::LogLevel::warn, "Could not open file {}", filename);
    }
    else {
        fprintf(file, "t0,dt_init,tmax,I0,total_pop,lambda,gamma\n");
        fprintf(file, "%.14f,%.14f,%.14f,%.14f,%.14f,%.14f,%.14f", t0, dt, tmax, I0, total_population, lambda, gamma);
        fclose(file);
    }
}
} // namespace setup
