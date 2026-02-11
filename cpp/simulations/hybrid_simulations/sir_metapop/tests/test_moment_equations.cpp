#include "smm_moments/simulation.h"
#include "smm_moments/model.h"
#include "smm_moments/parameters.h"
#include "smm_moments/closure_functions.h"
#include "simulations/hybrid_simulations/sir_metapop/library/moment_array.h"
#include "ode_sir/infection_state.h"
#include <cstddef>
#include <ostream>

namespace params
{
const double lambda1    = 0.0001;
const double lambda2    = 0.0002;
const double gamma      = 1. / 5.;
const double k_12_S     = 0.1;
const double k_21_S     = 0.025;
const double k_12_I     = 0.01;
const double k_21_I     = 0.003;
const double k_12_R     = 0.04;
const double k_21_R     = 0.0035;
const double total_pop1 = 10000;
const double total_pop2 = 10000;
const double I1         = 10;
const double I2         = 5;

} // namespace params

void moments_one_region(Eigen::Ref<const Eigen::VectorX<ScalarType>> y, ScalarType /*t*/,
                        Eigen::Ref<Eigen::VectorX<ScalarType>> dydt, double lambda, double gamma,
                        mio::smm_moments::Model<1, 4>& model)
{
    double mu_S = y[static_cast<size_t>(mio::osir::InfectionState::Susceptible)];
    double mu_I = y[static_cast<size_t>(mio::osir::InfectionState::Infected)];
    // current moment values
    // 2nd order
    double M_110 = y[model.moments.flatten_index({1, 1, 0}) + model.populations.get_num_compartments()];
    double M_101 = y[model.moments.flatten_index({1, 0, 1}) + model.populations.get_num_compartments()];
    double M_011 = y[model.moments.flatten_index({0, 1, 1}) + model.populations.get_num_compartments()];
    double M_200 = y[model.moments.flatten_index({2, 0, 0}) + model.populations.get_num_compartments()];
    double M_020 = y[model.moments.flatten_index({0, 2, 0}) + model.populations.get_num_compartments()];
    double M_002 = y[model.moments.flatten_index({0, 0, 2}) + model.populations.get_num_compartments()];
    // 3rd order
    double M_300 = y[model.moments.flatten_index({3, 0, 0}) + model.populations.get_num_compartments()];
    double M_030 = y[model.moments.flatten_index({0, 3, 0}) + model.populations.get_num_compartments()];
    double M_120 = y[model.moments.flatten_index({1, 2, 0}) + model.populations.get_num_compartments()];
    double M_102 = y[model.moments.flatten_index({1, 0, 2}) + model.populations.get_num_compartments()];
    double M_012 = y[model.moments.flatten_index({0, 1, 2}) + model.populations.get_num_compartments()];
    double M_111 = y[model.moments.flatten_index({1, 1, 1}) + model.populations.get_num_compartments()];
    double M_210 = y[model.moments.flatten_index({2, 1, 0}) + model.populations.get_num_compartments()];
    double M_201 = y[model.moments.flatten_index({2, 0, 1}) + model.populations.get_num_compartments()];
    double M_021 = y[model.moments.flatten_index({0, 2, 1}) + model.populations.get_num_compartments()];
    // equations
    dydt[static_cast<size_t>(mio::osir::InfectionState::Susceptible)] = -lambda * mu_S * mu_I - lambda * M_110;
    dydt[static_cast<size_t>(mio::osir::InfectionState::Infected)] =
        lambda * mu_S * mu_I - gamma * mu_I + lambda * M_110;
    dydt[static_cast<size_t>(mio::osir::InfectionState::Recovered)] = gamma * mu_I;
    // M110
    size_t index = model.moments.flatten_index({1, 1, 0}) + model.populations.get_num_compartments();
    dydt[index]  = -lambda * mu_S * mu_I + lambda * mu_I * M_200 - lambda * mu_S * M_020 +
                  (-lambda * mu_I + lambda * mu_S - gamma - lambda) * M_110 + lambda * M_210 - lambda * M_120;
    // M101
    index       = model.moments.flatten_index({1, 0, 1}) + model.populations.get_num_compartments();
    dydt[index] = -lambda * mu_I * M_101 - lambda * mu_S * M_011 + gamma * M_110 - lambda * M_111;
    // M011
    index = model.moments.flatten_index({0, 1, 1}) + model.populations.get_num_compartments();
    dydt[index] =
        -gamma * mu_I + lambda * mu_I * M_101 + gamma * M_020 + (lambda * mu_S - gamma) * M_011 + lambda * M_111;
    // M200
    index = model.moments.flatten_index({2, 0, 0}) + model.populations.get_num_compartments();
    dydt[index] =
        lambda * mu_S * mu_I - 2 * lambda * mu_I * M_200 + (lambda - 2 * lambda * mu_S) * M_110 - 2 * lambda * M_210;
    // M020
    index       = model.moments.flatten_index({0, 2, 0}) + model.populations.get_num_compartments();
    dydt[index] = lambda * mu_S * mu_I + gamma * mu_I + (2 * lambda * mu_S - 2 * gamma) * M_020 +
                  (2 * lambda * mu_I + lambda) * M_110 + 2 * lambda * M_120;
    // M002
    index       = model.moments.flatten_index({0, 0, 2}) + model.populations.get_num_compartments();
    dydt[index] = gamma * mu_I + 2 * gamma * M_011;
    // M300
    index       = model.moments.flatten_index({3, 0, 0}) + model.populations.get_num_compartments();
    dydt[index] = -lambda * mu_I * mu_S + (3 * lambda * mu_S - lambda) * M_110 +
                  (-3 * dydt[static_cast<size_t>(mio::osir::InfectionState::Susceptible)] - 3 * lambda * mu_I * mu_S +
                   3 * lambda * mu_I) *
                      M_200 +
                  (-3 * lambda * mu_S + 3 * lambda) * M_210 - 3 * lambda * mu_I * M_300;
    // M030
    index = model.moments.flatten_index({0, 3, 0}) + model.populations.get_num_compartments();
    dydt[index] =
        lambda * mu_I * mu_S - gamma * mu_I + (3 * lambda * mu_I * mu_S + 3 * lambda * mu_S) * M_020 +
        3 * lambda * mu_S * M_030 + (3 * lambda * mu_I + lambda) * M_110 + (3 * lambda * mu_I + 3 * lambda) * M_120 +
        (-3 * gamma * mu_I + 3 * gamma - 3 * dydt[static_cast<size_t>(mio::osir::InfectionState::Infected)]) * M_020 -
        3 * gamma * M_030;
    // M003
    index       = model.moments.flatten_index({0, 0, 3}) + model.populations.get_num_compartments();
    dydt[index] = gamma * mu_I +
                  (3 * gamma * mu_I - 3 * dydt[static_cast<size_t>(mio::osir::InfectionState::Recovered)]) * M_002 +
                  3 * gamma * M_011 + 3 * gamma * M_012;

    // M120
    index       = model.moments.flatten_index({1, 2, 0}) + model.populations.get_num_compartments();
    dydt[index] = -lambda * mu_I * mu_S +
                  (-lambda * mu_I * mu_S - 2 * lambda * mu_S -
                   dydt[static_cast<size_t>(mio::osir::InfectionState::Susceptible)]) *
                      M_020 -
                  lambda * mu_S * M_030 +
                  (2 * lambda * mu_I * mu_S - 2 * lambda * mu_I + lambda * mu_S - lambda - 2 * gamma * mu_I + gamma -
                   2 * dydt[static_cast<size_t>(mio::osir::InfectionState::Infected)]) *
                      M_110 +
                  (-lambda * mu_I + 2 * lambda * mu_S - 2 * lambda - 2 * gamma) * M_120 + lambda * mu_I * M_200 +
                  (2 * lambda * mu_I + lambda) * M_210;

    // M102
    index       = model.moments.flatten_index({1, 0, 2}) + model.populations.get_num_compartments();
    dydt[index] = (-lambda * mu_I * mu_S - dydt[static_cast<size_t>(mio::osir::InfectionState::Susceptible)]) * M_002 -
                  lambda * mu_S * M_012 - lambda * mu_I * M_102 +
                  (2 * gamma * mu_I - 2 * dydt[static_cast<size_t>(mio::osir::InfectionState::Recovered)]) * M_101 +
                  gamma * M_110 + 2 * gamma * M_111;

    // M012
    index = model.moments.flatten_index({0, 1, 2}) + model.populations.get_num_compartments();
    dydt[index] =
        -gamma * mu_I +
        (lambda * mu_I * mu_S - gamma * mu_I - dydt[static_cast<size_t>(mio::osir::InfectionState::Infected)]) * M_002 +
        (lambda * mu_S - gamma) * M_012 + lambda * mu_I * M_102 +
        (2 * gamma * mu_I - 2 * gamma - 2 * dydt[static_cast<size_t>(mio::osir::InfectionState::Recovered)]) * M_011 +
        gamma * M_020 + 2 * gamma * M_021;

    // M111
    index = model.moments.flatten_index({1, 1, 1}) + model.populations.get_num_compartments();
    dydt[index] =
        (-lambda * mu_I * mu_S - lambda * mu_S - dydt[static_cast<size_t>(mio::osir::InfectionState::Susceptible)]) *
            M_011 -
        lambda * mu_S * M_021 +
        (lambda * mu_I * mu_S - lambda * mu_I - gamma * mu_I -
         dydt[static_cast<size_t>(mio::osir::InfectionState::Infected)]) *
            M_101 +
        (-lambda * mu_I + lambda * mu_S - lambda - gamma) * M_111 + lambda * mu_I * M_201 +
        (gamma * mu_I - gamma - dydt[static_cast<size_t>(mio::osir::InfectionState::Recovered)]) * M_110 +
        gamma * M_120;

    // M210
    index       = model.moments.flatten_index({2, 1, 0}) + model.populations.get_num_compartments();
    dydt[index] = lambda * mu_I * mu_S + lambda * mu_S * M_020 +
                  (-2 * lambda * mu_I * mu_S + lambda * mu_I - 2 * lambda * mu_S + lambda -
                   2 * dydt[static_cast<size_t>(mio::osir::InfectionState::Susceptible)]) *
                      M_110 +
                  (-2 * lambda * mu_S + lambda) * M_120 +
                  (lambda * mu_I * mu_S - 2 * lambda * mu_I - gamma * mu_I -
                   dydt[static_cast<size_t>(mio::osir::InfectionState::Infected)]) *
                      M_200 +
                  (-2 * lambda * mu_I + lambda * mu_S - 2 * lambda - gamma) * M_210 + lambda * mu_I * M_300;

    // M201
    index       = model.moments.flatten_index({2, 0, 1}) + model.populations.get_num_compartments();
    dydt[index] = lambda * mu_S * M_011 +
                  (-2 * lambda * mu_I * mu_S + lambda * mu_I -
                   2 * dydt[static_cast<size_t>(mio::osir::InfectionState::Susceptible)]) *
                      M_101 +
                  (-2 * lambda * mu_S + lambda) * M_111 - 2 * lambda * mu_I * M_201 +
                  (gamma * mu_I - dydt[static_cast<size_t>(mio::osir::InfectionState::Recovered)]) * M_200 +
                  gamma * M_210;

    // M021
    index       = model.moments.flatten_index({0, 2, 1}) + model.populations.get_num_compartments();
    dydt[index] = gamma * mu_I +
                  (2 * lambda * mu_I * mu_S + lambda * mu_S - 2 * gamma * mu_I + gamma -
                   2 * dydt[static_cast<size_t>(mio::osir::InfectionState::Infected)]) *
                      M_011 +
                  (2 * lambda * mu_S - 2 * gamma) * M_021 + lambda * mu_I * M_101 +
                  (2 * lambda * mu_I + lambda) * M_111 +
                  (gamma * mu_I - 2 * gamma - dydt[static_cast<size_t>(mio::osir::InfectionState::Recovered)]) * M_020 +
                  gamma * M_030;
}

void test_one_region()
{
    std::cerr << "Running test one region..." << std::endl;
    // Initialize model and set parameters
    mio::smm_moments::Model<1, 4> model(&mio::smm_moments::truncation_closure<1>);
    model.populations[{mio::regions::Region(0), mio::osir::InfectionState::Susceptible}] =
        params::total_pop1 - params::I1;
    model.populations[{mio::regions::Region(0), mio::osir::InfectionState::Infected}]            = params::I1;
    model.populations[{mio::regions::Region(0), mio::osir::InfectionState::Recovered}]           = 0.0;
    model.parameters.template get<mio::smm_moments::TransmissionRate>()[mio::regions::Region(0)] = params::lambda1;
    model.parameters.template get<mio::smm_moments::RecoveryRate>()                              = params::gamma;

    Eigen::VectorX<ScalarType> y =
        Eigen::VectorX<ScalarType>::Zero(model.moments.moments().size() + model.populations.get_num_compartments());
    // mu_S, mu_I, mu_R
    y[0] = params::total_pop1 - params::I1;
    y[1] = params::I1;
    y[2] = 0.0;
    y[model.moments.flatten_index({0, 0, 0}) + model.populations.get_num_compartments()] = 1.0; // M000
    // second order moments
    y[model.moments.flatten_index({2, 0, 0}) + model.populations.get_num_compartments()] = 2.0; // M200
    y[model.moments.flatten_index({0, 2, 0}) + model.populations.get_num_compartments()] = 1.9; // M020
    y[model.moments.flatten_index({0, 0, 2}) + model.populations.get_num_compartments()] = 1.8; // M002
    y[model.moments.flatten_index({1, 1, 0}) + model.populations.get_num_compartments()] = 1.7; // M110
    y[model.moments.flatten_index({1, 0, 1}) + model.populations.get_num_compartments()] = 1.6; // M101
    y[model.moments.flatten_index({0, 1, 1}) + model.populations.get_num_compartments()] = 1.5; // M011
    // third order moments
    y[model.moments.flatten_index({3, 0, 0}) + model.populations.get_num_compartments()] = 1.4; // M300
    y[model.moments.flatten_index({0, 3, 0}) + model.populations.get_num_compartments()] = 1.3; // M030
    y[model.moments.flatten_index({0, 0, 3}) + model.populations.get_num_compartments()] = 1.2; // M003
    y[model.moments.flatten_index({1, 2, 0}) + model.populations.get_num_compartments()] = 1.1; // M120
    y[model.moments.flatten_index({1, 0, 2}) + model.populations.get_num_compartments()] = 1.0; // M102
    y[model.moments.flatten_index({0, 1, 2}) + model.populations.get_num_compartments()] = 0.9; // M012
    y[model.moments.flatten_index({1, 1, 1}) + model.populations.get_num_compartments()] = 0.8; // M111
    y[model.moments.flatten_index({2, 1, 0}) + model.populations.get_num_compartments()] = 0.7; // M210
    y[model.moments.flatten_index({2, 0, 1}) + model.populations.get_num_compartments()] = 0.6; // M201
    y[model.moments.flatten_index({0, 2, 1}) + model.populations.get_num_compartments()] = 0.5; // M021

    Eigen::VectorX<ScalarType> dydt1 = Eigen::VectorX<ScalarType>::Zero(y.size());
    Eigen::VectorX<ScalarType> dydt2 = Eigen::VectorX<ScalarType>::Zero(y.size());

    model.get_derivatives(y, 0.0, dydt1);
    moments_one_region(y, 0.0, dydt2, params::lambda1, params::gamma, model);

    double tol = 1e-10;
    for (size_t i = 0; i < static_cast<size_t>(dydt1.size()); ++i) {
        if (std::abs(dydt1[i] - dydt2[i]) > tol) {
            std::cerr << "Discrepancy found at index " << i << ": " << dydt1[i] << " vs " << dydt2[i] << std::endl;
            if (i >= 3) {
                auto multi_idx = model.moments.unflatten_index(i - model.populations.get_num_compartments());
                std::cerr << "  corresponding to moment M_";
                for (auto&& idx : multi_idx) {
                    std::cerr << idx;
                }
                std::cerr << std::endl;
            }
        }
    }
}

void moments_two_regions(Eigen::Ref<const Eigen::VectorX<ScalarType>> y, ScalarType /*t*/,
                         Eigen::Ref<Eigen::VectorX<ScalarType>> dydt, double lambda1, double lambda2, double gamma,
                         double k12_S, double k21_S, double k12_I, double k21_I, double k12_R, double k21_R,
                         mio::smm_moments::Model<2, 3>& model)
{
    double mu_S1 = y[static_cast<size_t>(mio::osir::InfectionState::Susceptible)];
    double mu_I1 = y[static_cast<size_t>(mio::osir::InfectionState::Infected)];
    double mu_R1 = y[static_cast<size_t>(mio::osir::InfectionState::Recovered)];
    double mu_S2 = y[static_cast<size_t>(mio::osir::InfectionState::Count) +
                     static_cast<size_t>(mio::osir::InfectionState::Susceptible)];
    double mu_I2 = y[static_cast<size_t>(mio::osir::InfectionState::Count) +
                     static_cast<size_t>(mio::osir::InfectionState::Infected)];
    double mu_R2 = y[static_cast<size_t>(mio::osir::InfectionState::Count) +
                     static_cast<size_t>(mio::osir::InfectionState::Recovered)];
    // current moment values
    // 2nd order
    double M_110000 = y[model.moments.flatten_index({1, 1, 0, 0, 0, 0}) + model.populations.get_num_compartments()];
    double M_101000 = y[model.moments.flatten_index({1, 0, 1, 0, 0, 0}) + model.populations.get_num_compartments()];
    double M_011000 = y[model.moments.flatten_index({0, 1, 1, 0, 0, 0}) + model.populations.get_num_compartments()];
    double M_200000 = y[model.moments.flatten_index({2, 0, 0, 0, 0, 0}) + model.populations.get_num_compartments()];
    double M_020000 = y[model.moments.flatten_index({0, 2, 0, 0, 0, 0}) + model.populations.get_num_compartments()];
    double M_002000 = y[model.moments.flatten_index({0, 0, 2, 0, 0, 0}) + model.populations.get_num_compartments()];
    double M_000110 = y[model.moments.flatten_index({0, 0, 0, 1, 1, 0}) + model.populations.get_num_compartments()];
    double M_000101 = y[model.moments.flatten_index({0, 0, 0, 1, 0, 1}) + model.populations.get_num_compartments()];
    double M_000011 = y[model.moments.flatten_index({0, 0, 0, 0, 1, 1}) + model.populations.get_num_compartments()];
    double M_000200 = y[model.moments.flatten_index({0, 0, 0, 2, 0, 0}) + model.populations.get_num_compartments()];
    double M_000020 = y[model.moments.flatten_index({0, 0, 0, 0, 2, 0}) + model.populations.get_num_compartments()];
    double M_000002 = y[model.moments.flatten_index({0, 0, 0, 0, 0, 2}) + model.populations.get_num_compartments()];
    double M_010100 = y[model.moments.flatten_index({0, 1, 0, 1, 0, 0}) + model.populations.get_num_compartments()];
    double M_100010 = y[model.moments.flatten_index({1, 0, 0, 0, 1, 0}) + model.populations.get_num_compartments()];
    double M_001100 = y[model.moments.flatten_index({0, 0, 1, 1, 0, 0}) + model.populations.get_num_compartments()];
    double M_100001 = y[model.moments.flatten_index({1, 0, 0, 0, 0, 1}) + model.populations.get_num_compartments()];
    double M_001010 = y[model.moments.flatten_index({0, 0, 1, 0, 1, 0}) + model.populations.get_num_compartments()];
    double M_010001 = y[model.moments.flatten_index({0, 1, 0, 0, 0, 1}) + model.populations.get_num_compartments()];
    double M_100100 = y[model.moments.flatten_index({1, 0, 0, 1, 0, 0}) + model.populations.get_num_compartments()];
    double M_010010 = y[model.moments.flatten_index({0, 1, 0, 0, 1, 0}) + model.populations.get_num_compartments()];
    double M_001001 = y[model.moments.flatten_index({0, 0, 1, 0, 0, 1}) + model.populations.get_num_compartments()];
    // 3rd order
    double M_210000 = y[model.moments.flatten_index({2, 1, 0, 0, 0, 0}) + model.populations.get_num_compartments()];
    double M_120000 = y[model.moments.flatten_index({1, 2, 0, 0, 0, 0}) + model.populations.get_num_compartments()];
    double M_111000 = y[model.moments.flatten_index({1, 1, 1, 0, 0, 0}) + model.populations.get_num_compartments()];
    double M_000210 = y[model.moments.flatten_index({0, 0, 0, 2, 1, 0}) + model.populations.get_num_compartments()];
    double M_000120 = y[model.moments.flatten_index({0, 0, 0, 1, 2, 0}) + model.populations.get_num_compartments()];
    double M_000111 = y[model.moments.flatten_index({0, 0, 0, 1, 1, 1}) + model.populations.get_num_compartments()];
    double M_110100 = y[model.moments.flatten_index({1, 1, 0, 1, 0, 0}) + model.populations.get_num_compartments()];
    double M_100110 = y[model.moments.flatten_index({1, 0, 0, 1, 1, 0}) + model.populations.get_num_compartments()];
    double M_110010 = y[model.moments.flatten_index({1, 1, 0, 0, 1, 0}) + model.populations.get_num_compartments()];
    double M_110001 = y[model.moments.flatten_index({1, 1, 0, 0, 0, 1}) + model.populations.get_num_compartments()];
    double M_010110 = y[model.moments.flatten_index({0, 1, 0, 1, 1, 0}) + model.populations.get_num_compartments()];
    double M_001110 = y[model.moments.flatten_index({0, 0, 1, 1, 1, 0}) + model.populations.get_num_compartments()];
    double M_011001 = y[model.moments.flatten_index({0, 1, 1, 0, 0, 1}) + model.populations.get_num_compartments()];
    //expected values
    dydt[static_cast<size_t>(mio::osir::InfectionState::Susceptible)] =
        -lambda1 * mu_S1 * mu_I1 - lambda1 * M_110000 - k12_S * mu_S1 + k21_S * mu_S2;
    dydt[static_cast<size_t>(mio::osir::InfectionState::Infected)] =
        lambda1 * mu_S1 * mu_I1 - gamma * mu_I1 + lambda1 * M_110000 - k12_I * mu_I1 + k21_I * mu_I2;
    dydt[static_cast<size_t>(mio::osir::InfectionState::Recovered)] = gamma * mu_I1 - k12_R * mu_R1 + k21_R * mu_R2;
    dydt[static_cast<size_t>(mio::osir::InfectionState::Count) +
         static_cast<size_t>(mio::osir::InfectionState::Susceptible)] =
        -lambda2 * mu_S2 * mu_I2 - lambda2 * M_000110 + k12_S * mu_S1 - k21_S * mu_S2;
    dydt[static_cast<size_t>(mio::osir::InfectionState::Count) +
         static_cast<size_t>(mio::osir::InfectionState::Infected)] =
        lambda2 * mu_S2 * mu_I2 - gamma * mu_I2 + lambda2 * M_000110 + k12_I * mu_I1 - k21_I * mu_I2;
    dydt[static_cast<size_t>(mio::osir::InfectionState::Count) +
         static_cast<size_t>(mio::osir::InfectionState::Recovered)] = gamma * mu_I2 + k12_R * mu_R1 - k21_R * mu_R2;
    // second-order moments
    //M110000
    size_t index = model.moments.flatten_index({1, 1, 0, 0, 0, 0}) + model.populations.get_num_compartments();
    dydt[index]  = -lambda1 * mu_S1 * mu_I1 + lambda1 * mu_I1 * M_200000 - lambda1 * mu_S1 * M_020000 +
                  (-lambda1 * mu_I1 + lambda1 * mu_S1 - gamma - lambda1 - k12_S - k12_I) * M_110000 +
                  lambda1 * M_210000 - lambda1 * M_120000 + k21_S * M_010100 + k21_I * M_100010;
    //M101000
    index       = model.moments.flatten_index({1, 0, 1, 0, 0, 0}) + model.populations.get_num_compartments();
    dydt[index] = gamma * M_110000 - lambda1 * mu_S1 * M_011000 + (-lambda1 * mu_I1 - k12_S - k12_R) * M_101000 -
                  lambda1 * M_111000 + k21_S * M_001100 + k21_R * M_100001;
    //M011000
    index       = model.moments.flatten_index({0, 1, 1, 0, 0, 0}) + model.populations.get_num_compartments();
    dydt[index] = -gamma * mu_I1 + lambda1 * mu_I1 * M_101000 + gamma * M_020000 +
                  (lambda1 * mu_S1 - gamma - k12_I - k12_R) * M_011000 + lambda1 * M_111000 + k21_I * M_001010 +
                  k21_R * M_010001;
    //M200000
    index       = model.moments.flatten_index({2, 0, 0, 0, 0, 0}) + model.populations.get_num_compartments();
    dydt[index] = lambda1 * mu_S1 * mu_I1 + k12_S * mu_S1 + k21_S * mu_S2 +
                  (-2 * lambda1 * mu_I1 - 2 * k12_S) * M_200000 + (lambda1 - 2 * lambda1 * mu_S1) * M_110000 -
                  2 * lambda1 * M_210000 + 2 * k21_S * M_100100;
    //M020000
    index       = model.moments.flatten_index({0, 2, 0, 0, 0, 0}) + model.populations.get_num_compartments();
    dydt[index] = lambda1 * mu_S1 * mu_I1 + gamma * mu_I1 + k12_I * mu_I1 + k21_I * mu_I2 +
                  (2 * lambda1 * mu_I1 + lambda1) * M_110000 +
                  (2 * lambda1 * mu_S1 - 2 * gamma - 2 * k12_I) * M_020000 + 2 * lambda1 * M_120000 +
                  2 * k21_I * M_010010;
    //M002000
    index       = model.moments.flatten_index({0, 0, 2, 0, 0, 0}) + model.populations.get_num_compartments();
    dydt[index] = gamma * mu_I1 + 2 * gamma * M_011000 + k12_R * mu_R1 + k21_R * mu_R2 - 2 * k12_R * M_002000 +
                  2 * k21_R * M_001001;
    //M000110
    index       = model.moments.flatten_index({0, 0, 0, 1, 1, 0}) + model.populations.get_num_compartments();
    dydt[index] = -lambda2 * mu_S2 * mu_I2 + lambda2 * mu_I2 * M_000200 +
                  (lambda2 * mu_S2 - lambda2 * mu_I2 - gamma - k21_S - k21_I - lambda2) * M_000110 -
                  lambda2 * mu_S2 * M_000020 + k12_S * M_100010 + k12_I * M_010100 + lambda2 * M_000210 -
                  lambda2 * M_000120;
    //M000101
    index       = model.moments.flatten_index({0, 0, 0, 1, 0, 1}) + model.populations.get_num_compartments();
    dydt[index] = (-lambda2 * mu_I2 - k21_R - k21_S) * M_000101 - lambda2 * mu_S2 * M_000011 + gamma * M_000110 +
                  k12_S * M_100001 + k12_R * M_001100 - lambda2 * M_000111;
    //M000011
    index       = model.moments.flatten_index({0, 0, 0, 0, 1, 1}) + model.populations.get_num_compartments();
    dydt[index] = -gamma * mu_I2 + lambda2 * mu_I2 * M_000101 + (lambda2 * mu_S2 - gamma - k21_I - k21_R) * M_000011 +
                  gamma * M_000020 + k12_I * M_010001 + k12_R * M_001010 + lambda2 * M_000111;
    //M000200
    index       = model.moments.flatten_index({0, 0, 0, 2, 0, 0}) + model.populations.get_num_compartments();
    dydt[index] = lambda2 * mu_S2 * mu_I2 + k12_S * mu_S1 + k21_S * mu_S2 +
                  (-2 * lambda2 * mu_I2 - 2 * k21_S) * M_000200 + (-2 * lambda2 * mu_S2 + lambda2) * M_000110 +
                  2 * k12_S * M_100100 - 2 * lambda2 * M_000210;
    //M000020
    index       = model.moments.flatten_index({0, 0, 0, 0, 2, 0}) + model.populations.get_num_compartments();
    dydt[index] = lambda2 * mu_S2 * mu_I2 + gamma * mu_I2 + k12_I * mu_I1 + k21_I * mu_I2 +
                  (2 * lambda2 * mu_I2 + lambda2) * M_000110 +
                  (2 * lambda2 * mu_S2 - 2 * gamma - 2 * k21_I) * M_000020 + 2 * k12_I * M_010010 +
                  2 * lambda2 * M_000120;
    //M000002
    index       = model.moments.flatten_index({0, 0, 0, 0, 0, 2}) + model.populations.get_num_compartments();
    dydt[index] = gamma * mu_I2 + k12_R * mu_R1 + k21_R * mu_R2 + 2 * gamma * M_000011 + 2 * k12_R * M_001001 -
                  2 * k21_R * M_000002;
    //M100100
    index       = model.moments.flatten_index({1, 0, 0, 1, 0, 0}) + model.populations.get_num_compartments();
    dydt[index] = -k12_S * mu_S1 - k21_S * mu_S2 + (-lambda1 * mu_I1 - k21_S - k12_S - lambda2 * mu_I2) * M_100100 -
                  lambda1 * mu_S1 * M_010100 - lambda2 * mu_S2 * M_100010 + k12_S * M_200000 + k21_S * M_000200 -
                  lambda1 * M_110100 - lambda2 * M_100110;
    //M100010
    index       = model.moments.flatten_index({1, 0, 0, 0, 1, 0}) + model.populations.get_num_compartments();
    dydt[index] = (-lambda1 * mu_I1 - k12_S - k21_I + lambda2 * mu_S2 - gamma) * M_100010 - lambda1 * mu_S1 * M_010010 +
                  k21_S * M_000110 + k12_I * M_110000 - lambda1 * M_110010 + lambda2 * mu_I2 * M_100100 +
                  lambda2 * M_100110;
    //M100001
    index       = model.moments.flatten_index({1, 0, 0, 0, 0, 1}) + model.populations.get_num_compartments();
    dydt[index] = (-lambda1 * mu_I1 - k12_S - k21_R) * M_100001 - lambda1 * mu_S1 * M_010001 + k21_S * M_000101 +
                  k12_R * M_101000 + gamma * M_100010 - lambda1 * M_110001;
    //M010100
    index       = model.moments.flatten_index({0, 1, 0, 1, 0, 0}) + model.populations.get_num_compartments();
    dydt[index] = lambda1 * mu_I1 * M_100100 + (lambda1 * mu_S1 - gamma - k21_S - k12_I - lambda2 * mu_I2) * M_010100 +
                  k12_S * M_110000 + k21_I * M_000110 - lambda2 * mu_S2 * M_010010 + lambda1 * M_110100 -
                  lambda2 * M_010110;
    //M010010
    index       = model.moments.flatten_index({0, 1, 0, 0, 1, 0}) + model.populations.get_num_compartments();
    dydt[index] = -k12_I * mu_I1 - k21_I * mu_I2 + lambda1 * mu_I1 * M_100010 + lambda2 * mu_I2 * M_010100 +
                  (lambda1 * mu_S1 + lambda2 * mu_S2 - 2. * gamma - k21_I - k12_I) * M_010010 + k12_I * M_020000 +
                  k21_I * M_000020 + lambda1 * M_110010 + lambda2 * M_010110;
    //M010001
    index       = model.moments.flatten_index({0, 1, 0, 0, 0, 1}) + model.populations.get_num_compartments();
    dydt[index] = lambda1 * mu_I1 * M_100001 + (lambda1 * mu_S1 - gamma - k12_I - k21_R) * M_010001 + gamma * M_010010 +
                  k21_I * M_000011 + k12_R * M_011000 + lambda1 * M_011001;
    //M001100
    index       = model.moments.flatten_index({0, 0, 1, 1, 0, 0}) + model.populations.get_num_compartments();
    dydt[index] = (-lambda2 * mu_I2 - k21_S - k12_R) * M_001100 - lambda2 * mu_S2 * M_001010 + gamma * M_010100 +
                  k12_S * M_101000 + k21_R * M_000101 - lambda2 * M_001110;
    //M001010
    index       = model.moments.flatten_index({0, 0, 1, 0, 1, 0}) + model.populations.get_num_compartments();
    dydt[index] = lambda2 * mu_I2 * M_001100 + (lambda2 * mu_S2 - gamma - k21_I - k12_R) * M_001010 + gamma * M_010010 +
                  k12_I * M_011000 + k21_R * M_000011 + lambda2 * M_001110;
    //M001001
    index       = model.moments.flatten_index({0, 0, 1, 0, 0, 1}) + model.populations.get_num_compartments();
    dydt[index] = -k12_R * mu_R1 - k21_R * mu_R2 + gamma * M_010001 + gamma * M_001010 + k12_R * M_002000 +
                  k21_R * M_000002 + (-k21_R - k12_R) * M_001001;
}

void test_two_regions()
{
    std::cerr << "Running test two regions..." << std::endl;
    // Initialize model and set parameters
    mio::smm_moments::Model<2, 3> model;
    // Region 0
    model.populations[{mio::regions::Region(0), mio::osir::InfectionState::Susceptible}] =
        params::total_pop1 - params::I1;
    model.populations[{mio::regions::Region(0), mio::osir::InfectionState::Infected}]              = params::I1;
    model.populations[{mio::regions::Region(0), mio::osir::InfectionState::Recovered}]             = 0.0;
    model.parameters.template get<mio::smm_moments::TransmissionRate>()[mio::regions::Region(0)]   = params::lambda1;
    model.parameters.template get<mio::smm_moments::RecoveryRate>()                                = params::gamma;
    model.parameters.template get<mio::smm_moments::TransitionRate>()[{
        mio::osir::InfectionState::Susceptible, mio::regions::Region(0), mio::regions::Region(1)}] = params::k_12_S;
    model.parameters.template get<mio::smm_moments::TransitionRate>()[{
        mio::osir::InfectionState::Infected, mio::regions::Region(0), mio::regions::Region(1)}]    = params::k_12_I;
    model.parameters.template get<mio::smm_moments::TransitionRate>()[{
        mio::osir::InfectionState::Recovered, mio::regions::Region(0), mio::regions::Region(1)}]   = params::k_12_R;
    //Region 1
    model.populations[{mio::regions::Region(1), mio::osir::InfectionState::Susceptible}] =
        params::total_pop2 - params::I2;
    model.populations[{mio::regions::Region(1), mio::osir::InfectionState::Infected}]              = params::I2;
    model.populations[{mio::regions::Region(1), mio::osir::InfectionState::Recovered}]             = 0.0;
    model.parameters.template get<mio::smm_moments::TransmissionRate>()[mio::regions::Region(1)]   = params::lambda2;
    model.parameters.template get<mio::smm_moments::TransitionRate>()[{
        mio::osir::InfectionState::Susceptible, mio::regions::Region(1), mio::regions::Region(0)}] = params::k_21_S;
    model.parameters.template get<mio::smm_moments::TransitionRate>()[{
        mio::osir::InfectionState::Infected, mio::regions::Region(1), mio::regions::Region(0)}]    = params::k_21_I;
    model.parameters.template get<mio::smm_moments::TransitionRate>()[{
        mio::osir::InfectionState::Recovered, mio::regions::Region(1), mio::regions::Region(0)}]   = params::k_21_R;

    Eigen::VectorX<ScalarType> y =
        Eigen::VectorX<ScalarType>::Zero(model.moments.moments().size() + model.populations.get_num_compartments());
    // mu_S, mu_I, mu_R
    // Region 0
    y[0] = params::total_pop1 - params::I1;
    y[1] = params::I1;
    y[2] = 0.0;
    // Region 1
    y[3] = params::total_pop2 - params::I2;
    y[4] = params::I2;
    y[5] = 0.0;
    y[model.moments.flatten_index({0, 0, 0, 0, 0, 0}) + model.populations.get_num_compartments()] = 1.0; // M000
    // second order moments
    // Region 0
    y[model.moments.flatten_index({2, 0, 0, 0, 0, 0}) + model.populations.get_num_compartments()] = 2.0; // M200000
    y[model.moments.flatten_index({0, 2, 0, 0, 0, 0}) + model.populations.get_num_compartments()] = 1.9; // M020000
    y[model.moments.flatten_index({0, 0, 2, 0, 0, 0}) + model.populations.get_num_compartments()] = 1.8; // M002000
    y[model.moments.flatten_index({1, 1, 0, 0, 0, 0}) + model.populations.get_num_compartments()] = 1.7; // M110000
    y[model.moments.flatten_index({1, 0, 1, 0, 0, 0}) + model.populations.get_num_compartments()] = 1.6; // M101000
    y[model.moments.flatten_index({0, 1, 1, 0, 0, 0}) + model.populations.get_num_compartments()] = 1.5; // M011000
    // Region 1
    y[model.moments.flatten_index({0, 0, 0, 2, 0, 0}) + model.populations.get_num_compartments()] = 6.0; // M000200
    y[model.moments.flatten_index({0, 0, 0, 0, 2, 0}) + model.populations.get_num_compartments()] = 5.9; // M000020
    y[model.moments.flatten_index({0, 0, 0, 0, 0, 2}) + model.populations.get_num_compartments()] = 5.8; // M000002
    y[model.moments.flatten_index({0, 0, 0, 1, 1, 0}) + model.populations.get_num_compartments()] = 5.7; // M000110
    y[model.moments.flatten_index({0, 0, 0, 1, 0, 1}) + model.populations.get_num_compartments()] = 5.6; // M000101
    y[model.moments.flatten_index({0, 0, 0, 0, 1, 1}) + model.populations.get_num_compartments()] = 5.5; // M000011
    // Mixed
    y[model.moments.flatten_index({1, 0, 0, 1, 0, 0}) + model.populations.get_num_compartments()] = -10.0; // M100100
    y[model.moments.flatten_index({1, 0, 0, 0, 1, 0}) + model.populations.get_num_compartments()] = -9.9; // M100010
    y[model.moments.flatten_index({1, 0, 0, 0, 0, 1}) + model.populations.get_num_compartments()] = -9.8; // M100001
    y[model.moments.flatten_index({0, 1, 0, 1, 0, 0}) + model.populations.get_num_compartments()] = -9.7; // M010100
    y[model.moments.flatten_index({0, 1, 0, 0, 1, 0}) + model.populations.get_num_compartments()] = -9.6; // M010010
    y[model.moments.flatten_index({0, 1, 0, 0, 0, 1}) + model.populations.get_num_compartments()] = -9.5; // M010001
    y[model.moments.flatten_index({0, 0, 1, 1, 0, 0}) + model.populations.get_num_compartments()] = -9.4; // M001100
    y[model.moments.flatten_index({0, 0, 1, 0, 1, 0}) + model.populations.get_num_compartments()] = -9.3; // M001010
    y[model.moments.flatten_index({0, 0, 1, 0, 0, 1}) + model.populations.get_num_compartments()] = -9.2; // M001001

    Eigen::VectorX<ScalarType> dydt1 = Eigen::VectorX<ScalarType>::Zero(y.size());
    Eigen::VectorX<ScalarType> dydt2 = Eigen::VectorX<ScalarType>::Zero(y.size());

    model.get_derivatives(y, 0.0, dydt1);
    moments_two_regions(y, 0.0, dydt2, params::lambda1, params::lambda2, params::gamma, params::k_12_S, params::k_21_S,
                        params::k_12_I, params::k_21_I, params::k_12_R, params::k_21_R, model);

    double tol = 1e-10;
    for (size_t i = 0; i < static_cast<size_t>(dydt1.size()); ++i) {
        if (std::abs(dydt1[i] - dydt2[i]) > tol) {
            std::cerr << "Discrepancy found at index " << i << ": " << dydt1[i] << " vs " << dydt2[i] << std::endl;
            if (i >= 6) {
                auto multi_idx = model.moments.unflatten_index(i - model.populations.get_num_compartments());
                std::cerr << "  corresponding to moment M_";
                for (auto&& idx : multi_idx) {
                    std::cerr << idx;
                }
                std::cerr << std::endl;
            }
        }
    }
}

int main()
{
    test_one_region();
    test_two_regions();
    return 0;
}
