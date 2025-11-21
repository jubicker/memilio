#include "hybrid_simulations/sir/moment_equation_array.h"
#include "memilio/math/euler.h"
#include "memilio/math/integrator.h"
#include "ode_sir/infection_state.h"

void second_and_third_order_moments(Eigen::Ref<const Eigen::VectorX<ScalarType>> y, ScalarType /*t*/,
                                    Eigen::Ref<Eigen::VectorX<ScalarType>> dydt, double lambda, double gamma,
                                    const MomentEquationArray<mio::osir::InfectionState, 1, 4>& values)
{
    double mu_S = y[static_cast<size_t>(mio::osir::InfectionState::Susceptible)];
    double mu_I = y[static_cast<size_t>(mio::osir::InfectionState::Infected)];
    // current moment values
    // 2nd order
    double M_110 = y[values.get_moment_index({MyIndex(1), MyIndex(1), MyIndex(0)})];
    double M_101 = y[values.get_moment_index({MyIndex(1), MyIndex(0), MyIndex(1)})];
    double M_011 = y[values.get_moment_index({MyIndex(0), MyIndex(1), MyIndex(1)})];
    double M_200 = y[values.get_moment_index({MyIndex(2), MyIndex(0), MyIndex(0)})];
    double M_020 = y[values.get_moment_index({MyIndex(0), MyIndex(2), MyIndex(0)})];
    double M_002 = y[values.get_moment_index({MyIndex(0), MyIndex(0), MyIndex(2)})];
    // 3rd order
    double M_300 = y[values.get_moment_index({MyIndex(3), MyIndex(0), MyIndex(0)})];
    double M_030 = y[values.get_moment_index({MyIndex(0), MyIndex(3), MyIndex(0)})];
    //double M_003 = y[values.get_moment_index({MyIndex(0), MyIndex(0), MyIndex(3)})];
    double M_120 = y[values.get_moment_index({MyIndex(1), MyIndex(2), MyIndex(0)})];
    double M_102 = y[values.get_moment_index({MyIndex(1), MyIndex(0), MyIndex(2)})];
    double M_012 = y[values.get_moment_index({MyIndex(0), MyIndex(1), MyIndex(2)})];
    double M_111 = y[values.get_moment_index({MyIndex(1), MyIndex(1), MyIndex(1)})];
    double M_210 = y[values.get_moment_index({MyIndex(2), MyIndex(1), MyIndex(0)})];
    double M_201 = y[values.get_moment_index({MyIndex(2), MyIndex(0), MyIndex(1)})];
    double M_021 = y[values.get_moment_index({MyIndex(0), MyIndex(2), MyIndex(1)})];
    // equations
    dydt[static_cast<size_t>(mio::osir::InfectionState::Susceptible)] = -lambda * mu_S * mu_I - lambda * M_110;
    dydt[static_cast<size_t>(mio::osir::InfectionState::Infected)] =
        lambda * mu_S * mu_I - gamma * mu_I + lambda * M_110;
    dydt[static_cast<size_t>(mio::osir::InfectionState::Recovered)] = gamma * mu_I;
    // M110
    size_t index = values.get_moment_index({MyIndex(1), MyIndex(1), MyIndex(0)});
    dydt[index]  = -lambda * mu_S * mu_I + lambda * mu_I * M_200 - lambda * mu_S * M_020 +
                  (-lambda * mu_I + lambda * mu_S - gamma - lambda) * M_110 + lambda * M_210 - lambda * M_120;
    // M101
    index       = values.get_moment_index({MyIndex(1), MyIndex(0), MyIndex(1)});
    dydt[index] = -lambda * mu_I * M_101 - lambda * mu_S * M_011 + gamma * M_110 - lambda * M_111;
    // M011
    index = values.get_moment_index({MyIndex(0), MyIndex(1), MyIndex(1)});
    dydt[index] =
        -gamma * mu_I + lambda * mu_I * M_101 + gamma * M_020 + (lambda * mu_S - gamma) * M_011 + lambda * M_111;
    // M200
    index = values.get_moment_index({MyIndex(2), MyIndex(0), MyIndex(0)});
    dydt[index] =
        lambda * mu_S * mu_I - 2 * lambda * mu_I * M_200 + (lambda - 2 * lambda * mu_S) * M_110 - 2 * lambda * M_210;
    // M020
    index       = values.get_moment_index({MyIndex(0), MyIndex(2), MyIndex(0)});
    dydt[index] = lambda * mu_S * mu_I + gamma * mu_I + (2 * lambda * mu_S - 2 * gamma) * M_020 +
                  (2 * lambda * mu_I + lambda) * M_110 + 2 * lambda * M_120;
    // M002
    index       = values.get_moment_index({MyIndex(0), MyIndex(0), MyIndex(2)});
    dydt[index] = gamma * mu_I + 2 * gamma * M_011;
    // M300
    index       = values.get_moment_index({MyIndex(3), MyIndex(0), MyIndex(0)});
    dydt[index] = -lambda * mu_I * mu_S + (3 * lambda * mu_S - lambda) * M_110 +
                  (-3 * dydt[static_cast<size_t>(mio::osir::InfectionState::Susceptible)] - 3 * lambda * mu_I * mu_S +
                   3 * lambda * mu_I) *
                      M_200 +
                  (-3 * lambda * mu_S + 3 * lambda) * M_210 - 3 * lambda * mu_I * M_300;
    // M030
    index = values.get_moment_index({MyIndex(0), MyIndex(3), MyIndex(0)});
    dydt[index] =
        lambda * mu_I * mu_S - gamma * mu_I + (3 * lambda * mu_I * mu_S + 3 * lambda * mu_S) * M_020 +
        3 * lambda * mu_S * M_030 + (3 * lambda * mu_I + lambda) * M_110 + (3 * lambda * mu_I + 3 * lambda) * M_120 +
        (-3 * gamma * mu_I + 3 * gamma - 3 * dydt[static_cast<size_t>(mio::osir::InfectionState::Infected)]) * M_020 -
        3 * gamma * M_030;
    // M003
    index       = values.get_moment_index({MyIndex(0), MyIndex(0), MyIndex(3)});
    dydt[index] = gamma * mu_I +
                  (3 * gamma * mu_I - 3 * dydt[static_cast<size_t>(mio::osir::InfectionState::Recovered)]) * M_002 +
                  3 * gamma * M_011 + 3 * gamma * M_012;

    // M120
    index       = values.get_moment_index({MyIndex(1), MyIndex(2), MyIndex(0)});
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
    index       = values.get_moment_index({MyIndex(1), MyIndex(0), MyIndex(2)});
    dydt[index] = (-lambda * mu_I * mu_S - dydt[static_cast<size_t>(mio::osir::InfectionState::Susceptible)]) * M_002 -
                  lambda * mu_S * M_012 - lambda * mu_I * M_102 +
                  (2 * gamma * mu_I - 2 * dydt[static_cast<size_t>(mio::osir::InfectionState::Recovered)]) * M_101 +
                  gamma * M_110 + 2 * gamma * M_111;

    // M012
    index = values.get_moment_index({MyIndex(0), MyIndex(1), MyIndex(2)});
    dydt[index] =
        -gamma * mu_I +
        (lambda * mu_I * mu_S - gamma * mu_I - dydt[static_cast<size_t>(mio::osir::InfectionState::Infected)]) * M_002 +
        (lambda * mu_S - gamma) * M_012 + lambda * mu_I * M_102 +
        (2 * gamma * mu_I - 2 * gamma - 2 * dydt[static_cast<size_t>(mio::osir::InfectionState::Recovered)]) * M_011 +
        gamma * M_020 + 2 * gamma * M_021;

    // M111
    index = values.get_moment_index({MyIndex(1), MyIndex(1), MyIndex(1)});
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
    index       = values.get_moment_index({MyIndex(2), MyIndex(1), MyIndex(0)});
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
    index       = values.get_moment_index({MyIndex(2), MyIndex(0), MyIndex(1)});
    dydt[index] = lambda * mu_S * M_011 +
                  (-2 * lambda * mu_I * mu_S + lambda * mu_I -
                   2 * dydt[static_cast<size_t>(mio::osir::InfectionState::Susceptible)]) *
                      M_101 +
                  (-2 * lambda * mu_S + lambda) * M_111 - 2 * lambda * mu_I * M_201 +
                  (gamma * mu_I - dydt[static_cast<size_t>(mio::osir::InfectionState::Recovered)]) * M_200 +
                  gamma * M_210;

    // M021
    index       = values.get_moment_index({MyIndex(0), MyIndex(2), MyIndex(1)});
    dydt[index] = gamma * mu_I +
                  (2 * lambda * mu_I * mu_S + lambda * mu_S - 2 * gamma * mu_I + gamma -
                   2 * dydt[static_cast<size_t>(mio::osir::InfectionState::Infected)]) *
                      M_011 +
                  (2 * lambda * mu_S - 2 * gamma) * M_021 + lambda * mu_I * M_101 +
                  (2 * lambda * mu_I + lambda) * M_111 +
                  (gamma * mu_I - 2 * gamma - dydt[static_cast<size_t>(mio::osir::InfectionState::Recovered)]) * M_020 +
                  gamma * M_030;
}

int main()
{
    const size_t closing_order = 4;
    const double total_pop     = 10000;
    const double I0            = 5;
    double lambda              = 0.0001;
    double gamma               = 1. / 5.;
    const double t0            = 0;
    const double dt            = 0.001;
    const double tmax          = 30;

    SIR::Moments<closing_order> moments1({total_pop - I0, I0, 0});
    mio::TimeSeries<ScalarType> result1(t0, moments1.values.get_values());

    SIR::Moments<closing_order> moments2({total_pop - I0, I0, 0});
    mio::TimeSeries<ScalarType> result2(t0, moments1.values.get_values());

    // First check whether get derivatives with closure = 3 and second order moments produce the same values
    Eigen::VectorX<ScalarType> y                                              = moments1.values.get_values();
    y[moments1.values.get_moment_index({MyIndex(1), MyIndex(1), MyIndex(0)})] = 0.1;
    y[moments1.values.get_moment_index({MyIndex(1), MyIndex(0), MyIndex(1)})] = 0.2;
    y[moments1.values.get_moment_index({MyIndex(0), MyIndex(1), MyIndex(1)})] = 0.3;
    y[moments1.values.get_moment_index({MyIndex(2), MyIndex(0), MyIndex(0)})] = 0.4;
    y[moments1.values.get_moment_index({MyIndex(0), MyIndex(2), MyIndex(0)})] = 0.5;
    y[moments1.values.get_moment_index({MyIndex(1), MyIndex(1), MyIndex(1)})] = -0.1;
    y[moments1.values.get_moment_index({MyIndex(1), MyIndex(0), MyIndex(2)})] = -0.2;
    y[moments1.values.get_moment_index({MyIndex(2), MyIndex(1), MyIndex(0)})] = 0.7;
    y[moments1.values.get_moment_index({MyIndex(3), MyIndex(0), MyIndex(0)})] = -0.4;
    y[moments1.values.get_moment_index({MyIndex(0), MyIndex(3), MyIndex(0)})] = 0.6;
    y[moments1.values.get_moment_index({MyIndex(1), MyIndex(2), MyIndex(0)})] = -0.6;
    Eigen::VectorX<ScalarType> dydt1 = Eigen::VectorX<ScalarType>::Zero(y.size());
    Eigen::VectorX<ScalarType> dydt2 = Eigen::VectorX<ScalarType>::Zero(y.size());
    second_and_third_order_moments(y, t0, dydt1, lambda, gamma, moments1.values);
    moments1.get_derivatives(y, t0, dydt2, closing_order, lambda, gamma);
    double tol = 1e-10;
    if (!dydt1.isApprox(dydt2, tol)) {
        std::cerr << "Error: Derivatives do not match for closure order 4!" << std::endl;
        std::cerr << "mu_S: " << dydt1[static_cast<size_t>(mio::osir::InfectionState::Susceptible)] << "  "
                  << dydt2[static_cast<size_t>(mio::osir::InfectionState::Susceptible)] << "\n";
        std::cerr << "mu_I: " << dydt1[static_cast<size_t>(mio::osir::InfectionState::Infected)] << "  "
                  << dydt2[static_cast<size_t>(mio::osir::InfectionState::Infected)] << "\n";
        std::cerr << "mu_R: " << dydt1[static_cast<size_t>(mio::osir::InfectionState::Recovered)] << "  "
                  << dydt2[static_cast<size_t>(mio::osir::InfectionState::Recovered)] << "\n";
        std::cerr << "M110: " << dydt1[moments1.values.get_moment_index({MyIndex(1), MyIndex(1), MyIndex(0)})] << "  "
                  << dydt2[moments1.values.get_moment_index({MyIndex(1), MyIndex(1), MyIndex(0)})] << "\n";
        std::cerr << "M101: " << dydt1[moments1.values.get_moment_index({MyIndex(1), MyIndex(0), MyIndex(1)})] << "  "
                  << dydt2[moments1.values.get_moment_index({MyIndex(1), MyIndex(0), MyIndex(1)})] << "\n";
        std::cerr << "M011: " << dydt1[moments1.values.get_moment_index({MyIndex(0), MyIndex(1), MyIndex(1)})] << "  "
                  << dydt2[moments1.values.get_moment_index({MyIndex(0), MyIndex(1), MyIndex(1)})] << "\n";
        std::cerr << "M200: " << dydt1[moments1.values.get_moment_index({MyIndex(2), MyIndex(0), MyIndex(0)})] << "  "
                  << dydt2[moments1.values.get_moment_index({MyIndex(2), MyIndex(0), MyIndex(0)})] << "\n";
        std::cerr << "M020: " << dydt1[moments1.values.get_moment_index({MyIndex(0), MyIndex(2), MyIndex(0)})] << "  "
                  << dydt2[moments1.values.get_moment_index({MyIndex(0), MyIndex(2), MyIndex(0)})] << "\n";
        std::cerr << "M002: " << dydt1[moments1.values.get_moment_index({MyIndex(0), MyIndex(0), MyIndex(2)})] << "  "
                  << dydt2[moments1.values.get_moment_index({MyIndex(0), MyIndex(0), MyIndex(2)})] << "\n";
        std::cerr << "M300: " << dydt1[moments1.values.get_moment_index({MyIndex(3), MyIndex(0), MyIndex(0)})] << "  "
                  << dydt2[moments1.values.get_moment_index({MyIndex(3), MyIndex(0), MyIndex(0)})] << "\n";
        std::cerr << "M030: " << dydt1[moments1.values.get_moment_index({MyIndex(0), MyIndex(3), MyIndex(0)})] << "  "
                  << dydt2[moments1.values.get_moment_index({MyIndex(0), MyIndex(3), MyIndex(0)})] << "\n";
        std::cerr << "M003: " << dydt1[moments1.values.get_moment_index({MyIndex(0), MyIndex(0), MyIndex(3)})] << "  "
                  << dydt2[moments1.values.get_moment_index({MyIndex(0), MyIndex(0), MyIndex(3)})] << "\n";
        std::cerr << "M120: " << dydt1[moments1.values.get_moment_index({MyIndex(1), MyIndex(2), MyIndex(0)})] << "  "
                  << dydt2[moments1.values.get_moment_index({MyIndex(1), MyIndex(2), MyIndex(0)})] << "\n";
        std::cerr << "M102: " << dydt1[moments1.values.get_moment_index({MyIndex(1), MyIndex(0), MyIndex(2)})] << "  "
                  << dydt2[moments1.values.get_moment_index({MyIndex(1), MyIndex(0), MyIndex(2)})] << "\n";
        std::cerr << "M012: " << dydt1[moments1.values.get_moment_index({MyIndex(0), MyIndex(1), MyIndex(2)})] << "  "
                  << dydt2[moments1.values.get_moment_index({MyIndex(0), MyIndex(1), MyIndex(2)})] << "\n";
        std::cerr << "M111: " << dydt1[moments1.values.get_moment_index({MyIndex(1), MyIndex(1), MyIndex(1)})] << "  "
                  << dydt2[moments1.values.get_moment_index({MyIndex(1), MyIndex(1), MyIndex(1)})] << "\n";
        std::cerr << "M210: " << dydt1[moments1.values.get_moment_index({MyIndex(2), MyIndex(1), MyIndex(0)})] << "  "
                  << dydt2[moments1.values.get_moment_index({MyIndex(2), MyIndex(1), MyIndex(0)})] << "\n";
        std::cerr << "M201: " << dydt1[moments1.values.get_moment_index({MyIndex(2), MyIndex(0), MyIndex(1)})] << "  "
                  << dydt2[moments1.values.get_moment_index({MyIndex(2), MyIndex(0), MyIndex(1)})] << "\n";
        std::cerr << "M021: " << dydt1[moments1.values.get_moment_index({MyIndex(0), MyIndex(2), MyIndex(1)})] << "  "
                  << dydt2[moments1.values.get_moment_index({MyIndex(0), MyIndex(2), MyIndex(1)})] << "\n";
        return 1;
    }
    else {
        std::cout << "Derivatives match for closure order 4." << std::endl;
    }

    mio::OdeIntegrator<ScalarType> integrator1(std::make_shared<mio::EulerIntegratorCore<ScalarType>>());

    double dt_integrator = dt;
    integrator1.advance(
        [closing_order, moments1, lambda, gamma](auto&& pop, auto&& t, auto&& dydt) {
            moments1.get_derivatives(pop, t, dydt, closing_order, lambda, gamma);
        },
        tmax, dt_integrator, result1);

    mio::OdeIntegrator<ScalarType> integrator2(std::make_shared<mio::EulerIntegratorCore<ScalarType>>());

    dt_integrator = dt;
    integrator2.advance(
        [moments2, lambda, gamma](auto&& pop, auto&& t, auto&& dydt) {
            second_and_third_order_moments(pop, t, dydt, lambda, gamma, moments2.values);
        },
        tmax, dt_integrator, result2);

    auto result_moments1 = moments1.get_moment_ts_and_names(result1, closing_order);
    auto result_moments2 = moments2.get_moment_ts_and_names(result2, closing_order);

    if (result_moments1.first.get_num_time_points() != result_moments2.first.get_num_time_points()) {
        std::cerr << "Error: Time series lengths do not match!" << std::endl;
        return 1;
    }
    if (!result_moments1.first.get_last_value().isApprox(result_moments2.first.get_last_value(), tol)) {
        std::cerr << "Error: Final values do not match!" << std::endl;
        std::cerr << "Result 1: " << result_moments1.first.get_last_value().transpose() << "\n";
        std::cerr << "Result 2: " << result_moments2.first.get_last_value().transpose() << "\n";
        return 1;
    }
    std::cout << "Integrating matches for closure order 4." << std::endl;
    return 0;
}
