#include "ode_sir/infection_state.h"
#include "smm/simulation.h"
#include "smm/model.h"
#include "smm/parameters.h"
#include "memilio/data/analyze_result.h"
#include "memilio/epidemiology/adoption_rate.h"
#include <string>

void run_smm_sim(int sim_num)
{
    using Model                 = mio::smm::Model<1, mio::osir::InfectionState>;
    ScalarType total_population = 1000;
    double num_I                = 1;

    Model model;
    model.populations[{mio::regions::Region(0), mio::osir::InfectionState::Infected}]    = num_I;
    model.populations[{mio::regions::Region(0), mio::osir::InfectionState::Recovered}]   = 0;
    model.populations[{mio::regions::Region(0), mio::osir::InfectionState::Susceptible}] = total_population - num_I;

    //Set lambda and tau
    double lambda = 0.002;
    double gamma  = 1. / 5.;

    // Second-order adoption rate is lambda * N
    std::vector<mio::AdoptionRate<mio::osir::InfectionState>> adoption_rates;
    adoption_rates.push_back({mio::osir::InfectionState::Susceptible,
                              mio::osir::InfectionState::Infected,
                              mio::regions::Region(0),
                              lambda * total_population,
                              {{mio::osir::InfectionState::Infected, 1.}}});
    adoption_rates.push_back({mio::osir::InfectionState::Infected,
                              mio::osir::InfectionState::Recovered,
                              mio::regions::Region(0),
                              gamma,
                              {}});

    model.parameters.get<mio::smm::AdoptionRates<mio::osir::InfectionState>>() = adoption_rates;

    ScalarType t0   = 0.;
    ScalarType tmax = 30.;
    ScalarType dt   = 0.1;

    auto sim = mio::smm::Simulation(model, t0, dt);
    sim.advance(tmax);

    std::string output_file =
        "/home/bick_ju/Documents/TemporalHybridModel/sim_outputs/SMM/" + std::to_string(sim_num) + "_comps.csv";

    int num_steps = static_cast<int>(tmax / dt) + 1;

    std::vector<double> interpolation_tps(num_steps);

    for (int i = 0; i < num_steps; ++i) {
        interpolation_tps[i] = i * dt;
    }
    auto result = mio::interpolate_simulation_result(sim.get_result(), interpolation_tps)
                      .export_csv(output_file, {"S", "I", "R"});
}

int main()
{
    int num_sims = 100;
    for (int sim = 0; sim < num_sims; ++sim) {
        run_smm_sim(sim);
    }
}
