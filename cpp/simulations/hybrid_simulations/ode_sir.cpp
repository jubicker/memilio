#include "memilio/config.h"
#include "memilio/data/analyze_result.h"
#include "ode_sir/infection_state.h"
#include "ode_sir/model.h"
#include "ode_sir/parameters.h"
#include "memilio/compartments/simulation.h"

int main()
{

    mio::set_log_level(mio::LogLevel::debug);

    ScalarType t0   = 0.;
    ScalarType tmax = 30.;
    ScalarType dt   = 0.1;

    ScalarType total_population = 1000;

    mio::log_info("Simulating SIR; t={} ... {} with dt = {}.", t0, tmax, dt);

    mio::osir::Model<ScalarType> model(1);

    model.populations[{mio::AgeGroup(0), mio::osir::InfectionState::Infected}]  = 1;
    model.populations[{mio::AgeGroup(0), mio::osir::InfectionState::Recovered}] = 0;
    model.populations[{mio::AgeGroup(0), mio::osir::InfectionState::Susceptible}] =
        total_population - model.populations[{mio::AgeGroup(0), mio::osir::InfectionState::Infected}] -
        model.populations[{mio::AgeGroup(0), mio::osir::InfectionState::Recovered}];
    //Set contact frequency to total population
    mio::ContactMatrixGroup& contact_matrix =
        model.parameters.get<mio::osir::ContactPatterns<ScalarType>>().get_cont_freq_mat();
    contact_matrix[0].get_baseline().setConstant(total_population);
    double lambda = 0.002;
    //rho = lambda*N
    model.parameters.set<mio::osir::TransmissionProbabilityOnContact<ScalarType>>(lambda);
    //gamma = 1 / t_I
    double gamma = 1. / 5.;
    model.parameters.set<mio::osir::TimeInfected<ScalarType>>(1. / gamma);

    model.check_constraints();

    auto sir = mio::simulate(t0, tmax, dt, model);

    std::string output_file = "/home/bick_ju/Documents/TemporalHybridModel/sim_outputs/ODE/sir.csv";

    int num_steps = static_cast<int>(tmax / dt) + 1;

    std::vector<double> interpolation_tps(num_steps);

    for (int i = 0; i < num_steps; ++i) {
        interpolation_tps[i] = i * dt;
    }
    auto result = mio::interpolate_simulation_result(sir, interpolation_tps).export_csv(output_file, {"S", "I", "R"});
}
