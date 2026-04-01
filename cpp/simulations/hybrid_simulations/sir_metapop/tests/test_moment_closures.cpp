#include "memilio/utils/logging.h"
#include "smm_moments/closure_functions.h"
#include "simulations/hybrid_simulations/sir_metapop/library/moment_array.h"
#include "ode_sir/infection_state.h"
#include "smm_moments/model.h"
#include <cmath>
#include <cstddef>
#include <numeric>
#include <ostream>

double get_E_XS_XI_XR(int closure_type, Eigen::Ref<const Eigen::VectorX<ScalarType>> y,
                      mio::smm_moments::Model<1, 3>& model)
{
    // current moment and expected values
    double mu_S = y[static_cast<size_t>(mio::osir::InfectionState::Susceptible)];
    double mu_I = y[static_cast<size_t>(mio::osir::InfectionState::Infected)];
    double mu_R = y[static_cast<size_t>(mio::osir::InfectionState::Recovered)];
    // 2nd order
    double M_110 = y[model.moments.flatten_index({1, 1, 0}) + model.populations.get_num_compartments()];
    double M_101 = y[model.moments.flatten_index({1, 0, 1}) + model.populations.get_num_compartments()];
    double M_011 = y[model.moments.flatten_index({0, 1, 1}) + model.populations.get_num_compartments()];
    double M_200 = y[model.moments.flatten_index({2, 0, 0}) + model.populations.get_num_compartments()];
    double M_020 = y[model.moments.flatten_index({0, 2, 0}) + model.populations.get_num_compartments()];
    double M_002 = y[model.moments.flatten_index({0, 0, 2}) + model.populations.get_num_compartments()];
    if (closure_type == 0) { // truncation
        return mu_S * mu_I * mu_R + mu_R * M_110 + mu_I * M_101 + mu_S * M_011;
    }
    else if (closure_type == 1) { // uncorrelation
        return mu_S * mu_I * mu_R;
    }
    else if (closure_type == 2) { // lognormal
        return std::exp(std::log((mu_S * mu_S) / (std::sqrt(M_200 + mu_S * mu_S))) +
                        std::log((mu_I * mu_I) / (std::sqrt(M_020 + mu_I * mu_I))) +
                        std::log((mu_R * mu_R) / (std::sqrt(M_002 + mu_R * mu_R))) +
                        0.5 * std::log((M_200 + mu_S * mu_S) / (mu_S * mu_S)) +
                        0.5 * std::log((M_020 + mu_I * mu_I) / (mu_I * mu_I)) +
                        0.5 * std::log((M_002 + mu_R * mu_R) / (mu_R * mu_R)) + std::log(1 + M_110 / (mu_S * mu_I)) +
                        std::log(1 + M_101 / (mu_S * mu_R)) + std::log(1 + M_011 / (mu_I * mu_R)));
    }
    else {
        mio::log_error("Unknown closure type {}.", closure_type);
    }
    return -1;
}

double get_E_XS3(int closure_type, Eigen::Ref<const Eigen::VectorX<ScalarType>> y, mio::smm_moments::Model<1, 3>& model)
{
    // current moment and expected values
    double mu_S = y[static_cast<size_t>(mio::osir::InfectionState::Susceptible)];
    // 2nd order
    double M_200 = y[model.moments.flatten_index({2, 0, 0}) + model.populations.get_num_compartments()];
    if (closure_type == 0) { // truncation
        return std::pow(mu_S, 3.0) + 3 * mu_S * M_200;
    }
    else if (closure_type == 1) { // uncorrelation
        return std::pow(mu_S, 3.0);
    }
    else if (closure_type == 2) { // lognormal
        return std::exp(3 * std::log((mu_S * mu_S) / (std::sqrt(M_200 + mu_S * mu_S))) +
                        4.5 * std::log((M_200 + mu_S * mu_S) / (mu_S * mu_S)));
    }
    else {
        mio::log_error("Unknown closure type {}.", closure_type);
    }
    return -1;
}

double get_E_XI3(int closure_type, Eigen::Ref<const Eigen::VectorX<ScalarType>> y, mio::smm_moments::Model<1, 3>& model)
{
    // current moment and expected values
    double mu_I = y[static_cast<size_t>(mio::osir::InfectionState::Infected)];
    // 2nd order
    double M_020 = y[model.moments.flatten_index({0, 2, 0}) + model.populations.get_num_compartments()];
    if (closure_type == 0) { // truncation
        return std::pow(mu_I, 3.0) + 3 * mu_I * M_020;
    }
    else if (closure_type == 1) { // uncorrelation
        return std::pow(mu_I, 3.0);
    }
    else if (closure_type == 2) { // lognormal
        return std::exp(3 * std::log((mu_I * mu_I) / (std::sqrt(M_020 + mu_I * mu_I))) +
                        4.5 * std::log((M_020 + mu_I * mu_I) / (mu_I * mu_I)));
    }
    else {
        mio::log_error("Unknown closure type {}.", closure_type);
    }
    return -1;
}

double get_E_XR3(int closure_type, Eigen::Ref<const Eigen::VectorX<ScalarType>> y, mio::smm_moments::Model<1, 3>& model)
{
    // current moment and expected values
    double mu_R = y[static_cast<size_t>(mio::osir::InfectionState::Recovered)];
    // 2nd order
    double M_002 = y[model.moments.flatten_index({0, 0, 2}) + model.populations.get_num_compartments()];
    if (closure_type == 0) { // truncation
        return std::pow(mu_R, 3.0) + 3 * mu_R * M_002;
    }
    else if (closure_type == 1) { // uncorrelation
        return std::pow(mu_R, 3.0);
    }
    else if (closure_type == 2) { // lognormal
        return std::exp(3 * std::log((mu_R * mu_R) / (std::sqrt(M_002 + mu_R * mu_R))) +
                        4.5 * std::log((M_002 + mu_R * mu_R) / (mu_R * mu_R)));
    }
    else {
        mio::log_error("Unknown closure type {}.", closure_type);
    }
    return -1;
}

double get_E_XS2_XI(int closure_type, Eigen::Ref<const Eigen::VectorX<ScalarType>> y,
                    mio::smm_moments::Model<1, 3>& model)
{
    // current moment and expected values
    double mu_S = y[static_cast<size_t>(mio::osir::InfectionState::Susceptible)];
    double mu_I = y[static_cast<size_t>(mio::osir::InfectionState::Infected)];
    // 2nd order
    double M_110 = y[model.moments.flatten_index({1, 1, 0}) + model.populations.get_num_compartments()];
    double M_200 = y[model.moments.flatten_index({2, 0, 0}) + model.populations.get_num_compartments()];
    double M_020 = y[model.moments.flatten_index({0, 2, 0}) + model.populations.get_num_compartments()];
    if (closure_type == 0) { // truncation
        return mu_S * mu_S * mu_I + 2 * mu_S * M_110 + mu_I * M_200;
    }
    else if (closure_type == 1) { // uncorrelation
        return mu_S * mu_S * mu_I;
    }
    else if (closure_type == 2) { // lognormal
        return std::exp(2 * std::log((mu_S * mu_S) / (std::sqrt(M_200 + mu_S * mu_S))) +
                        std::log((mu_I * mu_I) / (std::sqrt(M_020 + mu_I * mu_I))) +
                        2 * std::log((M_200 + mu_S * mu_S) / (mu_S * mu_S)) +
                        0.5 * std::log((M_020 + mu_I * mu_I) / (mu_I * mu_I)) +
                        2 * std::log(1 + M_110 / (mu_S * mu_I)));
    }
    else {
        mio::log_error("Unknown closure type {}.", closure_type);
    }
    return -1;
}

double get_E_XS2_XR(int closure_type, Eigen::Ref<const Eigen::VectorX<ScalarType>> y,
                    mio::smm_moments::Model<1, 3>& model)
{
    // current moment and expected values
    double mu_S = y[static_cast<size_t>(mio::osir::InfectionState::Susceptible)];
    double mu_R = y[static_cast<size_t>(mio::osir::InfectionState::Recovered)];
    // 2nd order
    double M_101 = y[model.moments.flatten_index({1, 0, 1}) + model.populations.get_num_compartments()];
    double M_200 = y[model.moments.flatten_index({2, 0, 0}) + model.populations.get_num_compartments()];
    double M_002 = y[model.moments.flatten_index({0, 0, 2}) + model.populations.get_num_compartments()];
    if (closure_type == 0) { // truncation
        return mu_S * mu_S * mu_R + 2 * mu_S * M_101 + mu_R * M_200;
    }
    else if (closure_type == 1) { // uncorrelation
        return mu_S * mu_S * mu_R;
    }
    else if (closure_type == 2) { // lognormal
        return std::exp(2 * std::log((mu_S * mu_S) / (std::sqrt(M_200 + mu_S * mu_S))) +
                        std::log((mu_R * mu_R) / (std::sqrt(M_002 + mu_R * mu_R))) +
                        2 * std::log((M_200 + mu_S * mu_S) / (mu_S * mu_S)) +
                        0.5 * std::log((M_002 + mu_R * mu_R) / (mu_R * mu_R)) +
                        2 * std::log(1 + M_101 / (mu_S * mu_R)));
    }
    else {
        mio::log_error("Unknown closure type {}.", closure_type);
    }
    return -1;
}

double get_E_XI2_XR(int closure_type, Eigen::Ref<const Eigen::VectorX<ScalarType>> y,
                    mio::smm_moments::Model<1, 3>& model)
{
    // current moment and expected values
    double mu_I = y[static_cast<size_t>(mio::osir::InfectionState::Infected)];
    double mu_R = y[static_cast<size_t>(mio::osir::InfectionState::Recovered)];
    // 2nd order
    double M_011 = y[model.moments.flatten_index({0, 1, 1}) + model.populations.get_num_compartments()];
    double M_020 = y[model.moments.flatten_index({0, 2, 0}) + model.populations.get_num_compartments()];
    double M_002 = y[model.moments.flatten_index({0, 0, 2}) + model.populations.get_num_compartments()];
    if (closure_type == 0) { // truncation
        return mu_I * mu_I * mu_R + 2 * mu_I * M_011 + mu_R * M_020;
    }
    else if (closure_type == 1) { // uncorrelation
        return mu_I * mu_I * mu_R;
    }
    else if (closure_type == 2) { // lognormal
        return std::exp(2 * std::log((mu_I * mu_I) / (std::sqrt(M_020 + mu_I * mu_I))) +
                        std::log((mu_R * mu_R) / (std::sqrt(M_002 + mu_R * mu_R))) +
                        2 * std::log((M_020 + mu_I * mu_I) / (mu_I * mu_I)) +
                        0.5 * std::log((M_002 + mu_R * mu_R) / (mu_R * mu_R)) +
                        2 * std::log(1 + M_011 / (mu_I * mu_R)));
    }
    else {
        mio::log_error("Unknown closure type {}.", closure_type);
    }
    return -1;
}

double get_E_XS_XI2(int closure_type, Eigen::Ref<const Eigen::VectorX<ScalarType>> y,
                    mio::smm_moments::Model<1, 3>& model)
{
    // current moment and expected values
    double mu_S = y[static_cast<size_t>(mio::osir::InfectionState::Susceptible)];
    double mu_I = y[static_cast<size_t>(mio::osir::InfectionState::Infected)];
    // 2nd order
    double M_110 = y[model.moments.flatten_index({1, 1, 0}) + model.populations.get_num_compartments()];
    double M_020 = y[model.moments.flatten_index({0, 2, 0}) + model.populations.get_num_compartments()];
    double M_200 = y[model.moments.flatten_index({2, 0, 0}) + model.populations.get_num_compartments()];
    if (closure_type == 0) { // truncation
        return mu_S * mu_I * mu_I + 2 * mu_I * M_110 + mu_S * M_020;
    }
    else if (closure_type == 1) { // uncorrelation
        return mu_S * mu_I * mu_I;
    }
    else if (closure_type == 2) { // lognormal
        return std::exp(2 * std::log((mu_I * mu_I) / (std::sqrt(M_020 + mu_I * mu_I))) +
                        std::log((mu_S * mu_S) / (std::sqrt(M_200 + mu_S * mu_S))) +
                        2 * std::log((M_020 + mu_I * mu_I) / (mu_I * mu_I)) +
                        0.5 * std::log((M_200 + mu_S * mu_S) / (mu_S * mu_S)) +
                        2 * std::log(1 + M_110 / (mu_S * mu_I)));
    }
    else {
        mio::log_error("Unknown closure type {}.", closure_type);
    }
    return -1;
}

double get_E_XS_XR2(int closure_type, Eigen::Ref<const Eigen::VectorX<ScalarType>> y,
                    mio::smm_moments::Model<1, 3>& model)
{
    // current moment and expected values
    double mu_S = y[static_cast<size_t>(mio::osir::InfectionState::Susceptible)];
    double mu_R = y[static_cast<size_t>(mio::osir::InfectionState::Recovered)];
    // 2nd order
    double M_101 = y[model.moments.flatten_index({1, 0, 1}) + model.populations.get_num_compartments()];
    double M_002 = y[model.moments.flatten_index({0, 0, 2}) + model.populations.get_num_compartments()];
    double M_200 = y[model.moments.flatten_index({2, 0, 0}) + model.populations.get_num_compartments()];
    if (closure_type == 0) { // truncation
        return mu_S * mu_R * mu_R + 2 * mu_R * M_101 + mu_S * M_002;
    }
    else if (closure_type == 1) { // uncorrelation
        return mu_S * mu_R * mu_R;
    }
    else if (closure_type == 2) { // lognormal
        return std::exp(2 * std::log((mu_R * mu_R) / (std::sqrt(M_002 + mu_R * mu_R))) +
                        std::log((mu_S * mu_S) / (std::sqrt(M_200 + mu_S * mu_S))) +
                        2 * std::log((M_002 + mu_R * mu_R) / (mu_R * mu_R)) +
                        0.5 * std::log((M_200 + mu_S * mu_S) / (mu_S * mu_S)) +
                        2 * std::log(1 + M_101 / (mu_S * mu_R)));
    }
    else {
        mio::log_error("Unknown closure type {}.", closure_type);
    }
    return -1;
}

double get_E_XI_XR2(int closure_type, Eigen::Ref<const Eigen::VectorX<ScalarType>> y,
                    mio::smm_moments::Model<1, 3>& model)
{
    // current moment and expected values
    double mu_I = y[static_cast<size_t>(mio::osir::InfectionState::Infected)];
    double mu_R = y[static_cast<size_t>(mio::osir::InfectionState::Recovered)];
    // 2nd order
    double M_011 = y[model.moments.flatten_index({0, 1, 1}) + model.populations.get_num_compartments()];
    double M_002 = y[model.moments.flatten_index({0, 0, 2}) + model.populations.get_num_compartments()];
    double M_020 = y[model.moments.flatten_index({0, 2, 0}) + model.populations.get_num_compartments()];
    if (closure_type == 0) { // truncation
        return mu_I * mu_R * mu_R + 2 * mu_R * M_011 + mu_I * M_002;
    }
    else if (closure_type == 1) { // uncorrelation
        return mu_I * mu_R * mu_R;
    }
    else if (closure_type == 2) { // lognormal
        return std::exp(2 * std::log((mu_R * mu_R) / (std::sqrt(M_002 + mu_R * mu_R))) +
                        std::log((mu_I * mu_I) / (std::sqrt(M_020 + mu_I * mu_I))) +
                        2 * std::log((M_002 + mu_R * mu_R) / (mu_R * mu_R)) +
                        0.5 * std::log((M_020 + mu_I * mu_I) / (mu_I * mu_I)) +
                        2 * std::log(1 + M_011 / (mu_I * mu_R)));
    }
    else {
        mio::log_error("Unknown closure type {}.", closure_type);
    }
    return -1;
}
using ClosureFunctionType =
    ScalarType (*)(std::array<int, static_cast<size_t>(mio::osir::InfectionState::Count) * 1> index,
                   Eigen::Ref<const Eigen::VectorX<ScalarType>> y,
                   const MomentArray<static_cast<size_t>(mio::osir::InfectionState::Count), 1, 3>& moments);

void test_closure(Eigen::Ref<const Eigen::VectorX<ScalarType>> y, mio::smm_moments::Model<1, 3>& model,
                  ClosureFunctionType closure_func, int closure_type)
{
    // current moment and expected values
    double mu_S = y[static_cast<size_t>(mio::osir::InfectionState::Susceptible)];
    double mu_I = y[static_cast<size_t>(mio::osir::InfectionState::Infected)];
    double mu_R = y[static_cast<size_t>(mio::osir::InfectionState::Recovered)];
    // 2nd order
    double M_110 = y[model.moments.flatten_index({1, 1, 0}) + model.populations.get_num_compartments()];
    double M_101 = y[model.moments.flatten_index({1, 0, 1}) + model.populations.get_num_compartments()];
    double M_011 = y[model.moments.flatten_index({0, 1, 1}) + model.populations.get_num_compartments()];
    double M_200 = y[model.moments.flatten_index({2, 0, 0}) + model.populations.get_num_compartments()];
    double M_020 = y[model.moments.flatten_index({0, 2, 0}) + model.populations.get_num_compartments()];
    double M_002 = y[model.moments.flatten_index({0, 0, 2}) + model.populations.get_num_compartments()];

    Eigen::VectorX<ScalarType> dydt1 = Eigen::VectorX<ScalarType>::Zero(y.size());
    Eigen::VectorX<ScalarType> dydt2 = Eigen::VectorX<ScalarType>::Zero(y.size());
    //M111
    size_t index = model.moments.flatten_index({1, 1, 1}) + model.populations.get_num_compartments();
    dydt1[index] =
        get_E_XS_XI_XR(closure_type, y, model) - mu_S * M_011 - mu_I * M_101 - mu_R * M_110 - mu_S * mu_I * mu_R;
    dydt2[index] = closure_func({1, 1, 1}, y, model.moments);
    //M300
    index        = model.moments.flatten_index({3, 0, 0}) + model.populations.get_num_compartments();
    dydt1[index] = get_E_XS3(closure_type, y, model) - std::pow(mu_S, 3.0) - 3 * mu_S * M_200;
    dydt2[index] = closure_func({3, 0, 0}, y, model.moments);
    //M030
    index        = model.moments.flatten_index({0, 3, 0}) + model.populations.get_num_compartments();
    dydt1[index] = get_E_XI3(closure_type, y, model) - std::pow(mu_I, 3.0) - 3 * mu_I * M_020;
    dydt2[index] = closure_func({0, 3, 0}, y, model.moments);
    //M003
    index        = model.moments.flatten_index({0, 0, 3}) + model.populations.get_num_compartments();
    dydt1[index] = get_E_XR3(closure_type, y, model) - std::pow(mu_R, 3.0) - 3 * mu_R * M_002;
    dydt2[index] = closure_func({0, 0, 3}, y, model.moments);
    //M210
    index        = model.moments.flatten_index({2, 1, 0}) + model.populations.get_num_compartments();
    dydt1[index] = get_E_XS2_XI(closure_type, y, model) - 2 * mu_S * M_110 - mu_I * mu_S * mu_S - mu_I * M_200;
    dydt2[index] = closure_func({2, 1, 0}, y, model.moments);
    //M201
    index        = model.moments.flatten_index({2, 0, 1}) + model.populations.get_num_compartments();
    dydt1[index] = get_E_XS2_XR(closure_type, y, model) - 2 * mu_S * M_101 - mu_R * mu_S * mu_S - mu_R * M_200;
    dydt2[index] = closure_func({2, 0, 1}, y, model.moments);
    //M021
    index        = model.moments.flatten_index({0, 2, 1}) + model.populations.get_num_compartments();
    dydt1[index] = get_E_XI2_XR(closure_type, y, model) - 2 * mu_I * M_011 - mu_R * mu_I * mu_I - mu_R * M_020;
    dydt2[index] = closure_func({0, 2, 1}, y, model.moments);
    //M120
    index        = model.moments.flatten_index({1, 2, 0}) + model.populations.get_num_compartments();
    dydt1[index] = get_E_XS_XI2(closure_type, y, model) - 2 * mu_I * M_110 - mu_S * mu_I * mu_I - mu_S * M_020;
    dydt2[index] = closure_func({1, 2, 0}, y, model.moments);
    //M102
    index        = model.moments.flatten_index({1, 0, 2}) + model.populations.get_num_compartments();
    dydt1[index] = get_E_XS_XR2(closure_type, y, model) - 2 * mu_R * M_101 - mu_S * mu_R * mu_R - mu_S * M_002;
    dydt2[index] = closure_func({1, 0, 2}, y, model.moments);
    //M012
    index        = model.moments.flatten_index({0, 1, 2}) + model.populations.get_num_compartments();
    dydt1[index] = get_E_XI_XR2(closure_type, y, model) - 2 * mu_R * M_011 - mu_I * mu_R * mu_R - mu_I * M_002;
    dydt2[index] = closure_func({0, 1, 2}, y, model.moments);

    double tol = 1e-10;
    for (size_t i = 0; i < static_cast<size_t>(dydt1.size()); ++i) {
        if (std::abs(dydt1[i] - dydt2[i]) > tol) {
            std::cerr << "Discrepancy found at index " << i << ": " << dydt1[i] << " vs " << dydt2[i] << std::endl;
            std::cerr << "Discrepancy is "
                      << ": " << dydt1[i] - dydt2[i] << std::endl;
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

void test_all_closures_3rd_order()
{
    // Initialize model and set parameters
    const size_t closure_order = 3;
    const size_t num_regions   = 1;
    mio::smm_moments::Model<num_regions, closure_order> model(
        &mio::smm_moments::truncation_closure<num_regions, closure_order>);
    model.populations[{mio::regions::Region(0), mio::osir::InfectionState::Susceptible}] = 10000 - 1;
    model.populations[{mio::regions::Region(0), mio::osir::InfectionState::Infected}]    = 1;
    model.populations[{mio::regions::Region(0), mio::osir::InfectionState::Recovered}]   = 0.0;

    Eigen::VectorX<ScalarType> y =
        Eigen::VectorX<ScalarType>::Zero(model.moments.moments().size() + model.populations.get_num_compartments());
    // mu_S, mu_I, mu_R
    y[0] = model.populations[{mio::regions::Region(0), mio::osir::InfectionState::Susceptible}];
    y[1] = model.populations[{mio::regions::Region(0), mio::osir::InfectionState::Infected}];
    y[2] = model.populations[{mio::regions::Region(0), mio::osir::InfectionState::Recovered}];
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

    std::cout << "Testing truncation closure...\n";
    test_closure(y, model, &mio::smm_moments::truncation_closure<num_regions, closure_order>, 0);
    std::cout << "Testing pairwise closure...\n";
    test_closure(y, model, &mio::smm_moments::pairapprox_closure<num_regions, closure_order>, 1);
    std::cout << "Testing lognormal closure...\n";
    test_closure(y, model, &mio::smm_moments::lognormal_closure<num_regions, closure_order>, 2);
    std::cout << "Testing lognormal zero inflation closure...\n";
    test_closure(y, model, &mio::smm_moments::lognormal_zero_inflation_closure<num_regions, closure_order>, 3);
}

double get_E_XS4(int closure_type, Eigen::Ref<const Eigen::VectorX<ScalarType>> y, mio::smm_moments::Model<1, 4>& model)
{
    // current moment and expected values
    double mu_S = y[static_cast<size_t>(mio::osir::InfectionState::Susceptible)];
    // 2nd order
    double M_200 = y[model.moments.flatten_index({2, 0, 0}) + model.populations.get_num_compartments()];
    // 3rd order
    double M_300 = y[model.moments.flatten_index({3, 0, 0}) + model.populations.get_num_compartments()];
    if (closure_type == 3) { // lognormal zero inflation
        double p = 1 - (std::pow(mu_S, 6.) * (1 + 3 * M_200) + std::pow(mu_S, 3.) * M_300) /
                           (std::pow(M_200 + mu_S * mu_S, 3.));
        double sigma_tilde_S =
            std::log((std::pow(mu_S, 4.) * (1 + 3 * M_200) + mu_S * M_300) / std::pow(M_200 + mu_S * mu_S, 2.));
        double mu_tilde_S = std::log(mu_S / (1 - p)) - 0.5 * sigma_tilde_S;
        return (1 - p) * std::exp(4 * mu_tilde_S + 8 * sigma_tilde_S);
    }
    else {
        mio::log_error("Unknown closure type {}.", closure_type);
    }
    return -1;
}

double get_E_XI4(int closure_type, Eigen::Ref<const Eigen::VectorX<ScalarType>> y, mio::smm_moments::Model<1, 4>& model)
{
    // current moment and expected values
    double mu_S = y[static_cast<size_t>(mio::osir::InfectionState::Susceptible)];
    double mu_I = y[static_cast<size_t>(mio::osir::InfectionState::Infected)];
    // 2nd order
    double M_200 = y[model.moments.flatten_index({2, 0, 0}) + model.populations.get_num_compartments()];
    double M_020 = y[model.moments.flatten_index({0, 2, 0}) + model.populations.get_num_compartments()];
    // 3rd order
    double M_300 = y[model.moments.flatten_index({3, 0, 0}) + model.populations.get_num_compartments()];
    double M_030 = y[model.moments.flatten_index({0, 3, 0}) + model.populations.get_num_compartments()];
    if (closure_type == 3) { // lognormal zero inflation
        double p = 1 - (std::pow(mu_S, 6.) * (1 + 3 * M_200) + std::pow(mu_S, 3.) * M_300) /
                           (std::pow(M_200 + mu_S * mu_S, 3.));
        double sigma_tilde_I =
            std::log((std::pow(mu_I, 4.) * (1 + 3 * M_020) + mu_I * M_030) / std::pow(M_020 + mu_I * mu_I, 2.));
        double mu_tilde_I = std::log(mu_I / (1 - p)) - 0.5 * sigma_tilde_I;
        return (1 - p) * std::exp(4 * mu_tilde_I + 8 * sigma_tilde_I);
    }
    else {
        mio::log_error("Unknown closure type {}.", closure_type);
    }
    return -1;
}

double get_E_XR4(int closure_type, Eigen::Ref<const Eigen::VectorX<ScalarType>> y, mio::smm_moments::Model<1, 4>& model)
{
    // current moment and expected values
    double mu_S = y[static_cast<size_t>(mio::osir::InfectionState::Susceptible)];
    double mu_R = y[static_cast<size_t>(mio::osir::InfectionState::Recovered)];
    // 2nd order
    double M_200 = y[model.moments.flatten_index({2, 0, 0}) + model.populations.get_num_compartments()];
    double M_002 = y[model.moments.flatten_index({0, 0, 2}) + model.populations.get_num_compartments()];
    // 3rd order
    double M_300 = y[model.moments.flatten_index({3, 0, 0}) + model.populations.get_num_compartments()];
    double M_003 = y[model.moments.flatten_index({0, 0, 3}) + model.populations.get_num_compartments()];
    if (closure_type == 3) { // lognormal zero inflation
        double p = 1 - (std::pow(mu_S, 6.) * (1 + 3 * M_200) + std::pow(mu_S, 3.) * M_300) /
                           (std::pow(M_200 + mu_S * mu_S, 3.));
        double sigma_tilde_R =
            std::log((std::pow(mu_R, 4.) * (1 + 3 * M_002) + mu_R * M_003) / std::pow(M_002 + mu_R * mu_R, 2.));
        double mu_tilde_R = std::log(mu_R / (1 - p)) - 0.5 * sigma_tilde_R;
        return (1 - p) * std::exp(4 * mu_tilde_R + 8 * sigma_tilde_R);
    }
    else {
        mio::log_error("Unknown closure type {}.", closure_type);
    }
    return -1;
}

double get_E_XS2XIXR(int closure_type, Eigen::Ref<const Eigen::VectorX<ScalarType>> y,
                     mio::smm_moments::Model<1, 4>& model)
{
    // current moment and expected values
    double mu_S = y[static_cast<size_t>(mio::osir::InfectionState::Susceptible)];
    double mu_I = y[static_cast<size_t>(mio::osir::InfectionState::Infected)];
    double mu_R = y[static_cast<size_t>(mio::osir::InfectionState::Recovered)];
    // 2nd order
    double M_200 = y[model.moments.flatten_index({2, 0, 0}) + model.populations.get_num_compartments()];
    double M_020 = y[model.moments.flatten_index({0, 2, 0}) + model.populations.get_num_compartments()];
    double M_002 = y[model.moments.flatten_index({0, 0, 2}) + model.populations.get_num_compartments()];
    double M_110 = y[model.moments.flatten_index({1, 1, 0}) + model.populations.get_num_compartments()];
    double M_101 = y[model.moments.flatten_index({1, 0, 1}) + model.populations.get_num_compartments()];
    double M_011 = y[model.moments.flatten_index({0, 1, 1}) + model.populations.get_num_compartments()];
    // 3rd order
    double M_300 = y[model.moments.flatten_index({3, 0, 0}) + model.populations.get_num_compartments()];
    double M_030 = y[model.moments.flatten_index({0, 3, 0}) + model.populations.get_num_compartments()];
    double M_003 = y[model.moments.flatten_index({0, 0, 3}) + model.populations.get_num_compartments()];
    if (closure_type == 3) { // lognormal zero inflation
        double p = 1 - (std::pow(mu_S, 6.) * (1 + 3 * M_200) + std::pow(mu_S, 3.) * M_300) /
                           (std::pow(M_200 + mu_S * mu_S, 3.));
        double sigma_tilde_S =
            std::log((std::pow(mu_S, 4.) * (1 + 3 * M_200) + mu_S * M_300) / std::pow(M_200 + mu_S * mu_S, 2.));
        double mu_tilde_S = std::log(mu_S / (1 - p)) - 0.5 * sigma_tilde_S;
        double sigma_tilde_I =
            std::log((std::pow(mu_I, 4.) * (1 + 3 * M_020) + mu_I * M_030) / std::pow(M_020 + mu_I * mu_I, 2.));
        double mu_tilde_I = std::log(mu_I / (1 - p)) - 0.5 * sigma_tilde_I;
        double sigma_tilde_R =
            std::log((std::pow(mu_R, 4.) * (1 + 3 * M_002) + mu_R * M_003) / std::pow(M_002 + mu_R * mu_R, 2.));
        double mu_tilde_R = std::log(mu_R / (1 - p)) - 0.5 * sigma_tilde_R;
        double c_SI       = std::log((mu_S * mu_I + M_110) / (1 - p)) - mu_tilde_S - 0.5 * sigma_tilde_S - mu_tilde_I -
                      0.5 * sigma_tilde_I;
        double c_SR = std::log((mu_S * mu_R + M_101) / (1 - p)) - mu_tilde_S - 0.5 * sigma_tilde_S - mu_tilde_R -
                      0.5 * sigma_tilde_R;
        double c_IR = std::log((mu_I * mu_R + M_011) / (1 - p)) - mu_tilde_I - 0.5 * sigma_tilde_I - mu_tilde_R -
                      0.5 * sigma_tilde_R;
        return (1 - p) * std::exp(2 * mu_tilde_S + 2 * sigma_tilde_S + mu_tilde_I + 0.5 * sigma_tilde_I + mu_tilde_R +
                                  0.5 * sigma_tilde_R + 2 * c_SI + 2 * c_SR + c_IR);
    }
    else {
        mio::log_error("Unknown closure type {}.", closure_type);
    }
    return -1;
}

double get_E_XSXI3(int closure_type, Eigen::Ref<const Eigen::VectorX<ScalarType>> y,
                   mio::smm_moments::Model<1, 4>& model)
{
    // current moment and expected values
    double mu_S = y[static_cast<size_t>(mio::osir::InfectionState::Susceptible)];
    double mu_I = y[static_cast<size_t>(mio::osir::InfectionState::Infected)];
    // 2nd order
    double M_200 = y[model.moments.flatten_index({2, 0, 0}) + model.populations.get_num_compartments()];
    double M_020 = y[model.moments.flatten_index({0, 2, 0}) + model.populations.get_num_compartments()];
    double M_110 = y[model.moments.flatten_index({1, 1, 0}) + model.populations.get_num_compartments()];
    // 3rd order
    double M_300 = y[model.moments.flatten_index({3, 0, 0}) + model.populations.get_num_compartments()];
    double M_030 = y[model.moments.flatten_index({0, 3, 0}) + model.populations.get_num_compartments()];
    if (closure_type == 3) { // lognormal zero inflation
        double p = 1 - (std::pow(mu_S, 6.) * (1 + 3 * M_200) + std::pow(mu_S, 3.) * M_300) /
                           (std::pow(M_200 + mu_S * mu_S, 3.));
        double sigma_tilde_S =
            std::log((std::pow(mu_S, 4.) * (1 + 3 * M_200) + mu_S * M_300) / std::pow(M_200 + mu_S * mu_S, 2.));
        double mu_tilde_S = std::log(mu_S / (1 - p)) - 0.5 * sigma_tilde_S;
        double sigma_tilde_I =
            std::log((std::pow(mu_I, 4.) * (1 + 3 * M_020) + mu_I * M_030) / std::pow(M_020 + mu_I * mu_I, 2.));
        double mu_tilde_I = std::log(mu_I / (1 - p)) - 0.5 * sigma_tilde_I;

        double c_SI = std::log((mu_S * mu_I + M_110) / (1 - p)) - mu_tilde_S - 0.5 * sigma_tilde_S - mu_tilde_I -
                      0.5 * sigma_tilde_I;
        return (1 - p) * std::exp(mu_tilde_S + 0.5 * sigma_tilde_S + 3 * mu_tilde_I + 4.5 * sigma_tilde_I + 3 * c_SI);
    }
    else {
        mio::log_error("Unknown closure type {}.", closure_type);
    }
    return -1;
}
using ClosureFunctionTypeOrder4 =
    ScalarType (*)(std::array<int, static_cast<size_t>(mio::osir::InfectionState::Count) * 1> index,
                   Eigen::Ref<const Eigen::VectorX<ScalarType>> y,
                   const MomentArray<static_cast<size_t>(mio::osir::InfectionState::Count), 1, 4>& moments);

void test_closure(Eigen::Ref<const Eigen::VectorX<ScalarType>> y, mio::smm_moments::Model<1, 4>& model,
                  ClosureFunctionTypeOrder4 closure_func, int closure_type)
{
    // current moment and expected values
    double mu_S = y[static_cast<size_t>(mio::osir::InfectionState::Susceptible)];
    double mu_I = y[static_cast<size_t>(mio::osir::InfectionState::Infected)];
    double mu_R = y[static_cast<size_t>(mio::osir::InfectionState::Recovered)];
    // 2nd order
    double M_200 = y[model.moments.flatten_index({2, 0, 0}) + model.populations.get_num_compartments()];
    double M_020 = y[model.moments.flatten_index({0, 2, 0}) + model.populations.get_num_compartments()];
    double M_002 = y[model.moments.flatten_index({0, 0, 2}) + model.populations.get_num_compartments()];
    double M_110 = y[model.moments.flatten_index({1, 1, 0}) + model.populations.get_num_compartments()];
    double M_101 = y[model.moments.flatten_index({1, 0, 1}) + model.populations.get_num_compartments()];
    double M_011 = y[model.moments.flatten_index({0, 1, 1}) + model.populations.get_num_compartments()];
    // 3rd order
    double M_300 = y[model.moments.flatten_index({3, 0, 0}) + model.populations.get_num_compartments()];
    double M_030 = y[model.moments.flatten_index({0, 3, 0}) + model.populations.get_num_compartments()];
    double M_003 = y[model.moments.flatten_index({0, 0, 3}) + model.populations.get_num_compartments()];
    double M_210 = y[model.moments.flatten_index({2, 1, 0}) + model.populations.get_num_compartments()];
    double M_201 = y[model.moments.flatten_index({2, 0, 1}) + model.populations.get_num_compartments()];
    double M_120 = y[model.moments.flatten_index({1, 2, 0}) + model.populations.get_num_compartments()];
    double M_111 = y[model.moments.flatten_index({1, 1, 1}) + model.populations.get_num_compartments()];

    Eigen::VectorX<ScalarType> dydt1 = Eigen::VectorX<ScalarType>::Zero(y.size());
    Eigen::VectorX<ScalarType> dydt2 = Eigen::VectorX<ScalarType>::Zero(y.size());
    //M400
    size_t index = model.moments.flatten_index({4, 0, 0}) + model.populations.get_num_compartments();
    dydt1[index] = get_E_XS4(closure_type, y, model) - std::pow(mu_S, 4.0) - 6 * mu_S * mu_S * M_200 - 4 * mu_S * M_300;
    dydt2[index] = closure_func({4, 0, 0}, y, model.moments);
    //M040
    index        = model.moments.flatten_index({0, 4, 0}) + model.populations.get_num_compartments();
    dydt1[index] = get_E_XI4(closure_type, y, model) - std::pow(mu_I, 4.0) - 6 * mu_I * mu_I * M_020 - 4 * mu_I * M_030;
    dydt2[index] = closure_func({0, 4, 0}, y, model.moments);
    //M004
    index        = model.moments.flatten_index({0, 0, 4}) + model.populations.get_num_compartments();
    dydt1[index] = get_E_XR4(closure_type, y, model) - std::pow(mu_R, 4.0) - 6 * mu_R * mu_R * M_002 - 4 * mu_R * M_003;
    dydt2[index] = closure_func({0, 0, 4}, y, model.moments);
    //M211
    index        = model.moments.flatten_index({2, 1, 1}) + model.populations.get_num_compartments();
    dydt1[index] = get_E_XS2XIXR(closure_type, y, model) - mu_S * mu_S * M_011 - 2 * mu_S * M_111 -
                   2 * mu_S * mu_I * M_101 - mu_I * M_201 - mu_S * mu_S * mu_I * mu_R - 2 * mu_S * mu_R * M_110 -
                   mu_I * mu_R * M_200 - mu_R * M_210;
    dydt2[index] = closure_func({2, 1, 1}, y, model.moments);
    //M130
    index        = model.moments.flatten_index({1, 3, 0}) + model.populations.get_num_compartments();
    dydt1[index] = get_E_XSXI3(closure_type, y, model) - mu_S * std::pow(mu_I, 3.) - mu_S * M_030 -
                   3 * mu_I * mu_I * M_110 - 3 * mu_S * mu_I * M_020 - 3 * mu_I * M_120;
    dydt2[index] = closure_func({1, 3, 0}, y, model.moments);

    double tol = 1e-10;
    for (size_t i = 0; i < static_cast<size_t>(dydt1.size()); ++i) {
        if (std::abs(dydt1[i] - dydt2[i]) > tol) {
            std::cerr << "Discrepancy found at index " << i << ": " << dydt1[i] << " vs " << dydt2[i] << std::endl;
            std::cerr << "Discrepancy is "
                      << ": " << dydt1[i] - dydt2[i] << std::endl;
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

void test_closure_4th_order()
{
    // Initialize model and set parameters
    const size_t closure_order = 4;
    const size_t num_regions   = 1;
    mio::smm_moments::Model<num_regions, closure_order> model(
        &mio::smm_moments::truncation_closure<num_regions, closure_order>);
    model.populations[{mio::regions::Region(0), mio::osir::InfectionState::Susceptible}] = 10000 - 2;
    model.populations[{mio::regions::Region(0), mio::osir::InfectionState::Infected}]    = 1;
    model.populations[{mio::regions::Region(0), mio::osir::InfectionState::Recovered}]   = 1;

    Eigen::VectorX<ScalarType> y =
        Eigen::VectorX<ScalarType>::Zero(model.moments.moments().size() + model.populations.get_num_compartments());
    // mu_S, mu_I, mu_R
    y[0] = model.populations[{mio::regions::Region(0), mio::osir::InfectionState::Susceptible}];
    y[1] = model.populations[{mio::regions::Region(0), mio::osir::InfectionState::Infected}];
    y[2] = model.populations[{mio::regions::Region(0), mio::osir::InfectionState::Recovered}];
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

    std::cout << "Testing lognormal zero inflation closure...\n";
    test_closure(y, model, &mio::smm_moments::lognormal_zero_inflation_closure<num_regions, closure_order>, 3);
}

int main()
{
    //test_all_closures_3rd_order();
    test_closure_4th_order();

    return 0;
}
