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

#ifndef MOMENTS_CLOSURE_FUNCTIONS_H
#define MOMENTS_CLOSURE_FUNCTIONS_H

#include "memilio/config.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/logging.h"
#include "models/ode_sir/infection_state.h"
#include "memilio/math/eigen.h"
#include "simulations/hybrid_simulations/sir_metapop/library/moment_array.h"
#include <array>
#include <boost/math/special_functions/math_fwd.hpp>
#include <cmath>
#include <cstddef>
#include <numeric>

namespace mio
{
namespace smm_moments
{

/**
 * This file contains various closure functions. The functions provide approximations for the moment given by index.
 */

/**
 * @brief Implements a truncation closure i.e. the given moment is set to zero.
 * @tparam NumRegions Number of regions.
 * @tparam ClosureOrder Order to be closed.
 * @param[in] index Index of the central moment that should be approximated.
 * @param[in] y Not used.
 * @param[in] moments Not used.
 */
template <size_t NumRegions, size_t ClosureOrder>
ScalarType truncation_closure(
    std::array<int, static_cast<size_t>(osir::InfectionState::Count) * NumRegions> index,
    Eigen::Ref<const Eigen::VectorX<ScalarType>> y,
    const MomentArray<static_cast<size_t>(osir::InfectionState::Count), NumRegions, ClosureOrder>& moments)
{
    int order = std::accumulate(index.begin(), index.end(), 0);
    if (order != ClosureOrder) {
        mio::log_error("The order of the moment that should be approximated is unequal to the given closure order. "
                       "Moment order is {} and closure order is {}.",
                       order, ClosureOrder);
    }
    mio::unused(index, y, moments);
    return 0;
}

/**
 * @brief Returns the first parameter of a lognormally distributed random variable given its mean and variance.
 * @param[in] var Variance of the lognormally distributed random variable.
 * @param[in] mean Mean of the lognormally distributed random variable.
 */
ScalarType get_log_mu(ScalarType var, ScalarType mean);

/**
 * @brief Returns the second parameter of a lognormally distributed random variable given its mean and variance.
 * @param[in] var Variance of the lognormally distributed random variable.
 * @param[in] mean Mean of the lognormally distributed random variable.
 */
ScalarType get_log_sigma(ScalarType var, ScalarType mean);

/**
 * @brief Returns Cov(ln(X), ln(Y)) of two lognormally distributed random variables X and Y given their covariance and means.
 * @param[in] cov Covariance of the lognormally distributed random variables.
 * @param[in] mean1 Mean of the first lognormally distributed random variable.
 * @param[in] mean2 Mean of the second lognormally distributed random variable.
 */
ScalarType get_log_cov(ScalarType cov, ScalarType mean1, ScalarType mean2);

/**
 * @brief Calculates raw moment E[X_S^lS * X_I^lI * X_R^lR] for one region from central moments.
 * @param[in] lS Power of X_S in the raw moment.
 * @param[in] lI Power of X_I in the raw moment.
 * @param[in] lR Power of X_R in the raw moment.
 * @param[in] y Current value of expected value and all central moments.
 * @param[in] moments Moment array. Is only used to get the correct flat index of a moment in y.
 */
template <size_t NumRegions, size_t ClosureOrder>
ScalarType get_E_XS_lS_XI_lI_XR_lR(
    size_t lS, size_t lI, size_t lR, Eigen::Ref<const Eigen::VectorX<ScalarType>> y,
    const MomentArray<static_cast<size_t>(osir::InfectionState::Count), NumRegions, ClosureOrder>& moments)
{
    ScalarType res = 0;
    for (int iS = 0; iS <= (int)lS; ++iS) {
        for (int iI = 0; iI <= (int)lI; ++iI) {
            for (int iR = 0; iR <= (int)lR; ++iR) {
                res += boost::math::binomial_coefficient<double>(lS, iS) *
                       boost::math::binomial_coefficient<double>(lI, iI) *
                       boost::math::binomial_coefficient<double>(lR, iR) * std::pow(y[0], lS - iS) *
                       std::pow(y[1], lI - iI) * std::pow(y[2], lR - iR) *
                       y[moments.flatten_index({iS, iI, iR}) + static_cast<size_t>(osir::InfectionState::Count)];
            }
        }
    }
    return res;
}

/**
 * @brief Calculates central moment given by index for one region from raw moments.
 * @tparam NumRegions Number of regions.
 * @tparam ClosureOrder Order to be closed.
 * @param[in] index Index of the moment that should be approximated.
 * @param[in] y Current value of expected value and all central moments.
 * @param[in] moments Moment array. Is only used to get the correct flat index of a moment in y.
 * @param[in] E_XS_rS_XI_rI_XR_rR Approximation of raw moment corresponding to central moment that should be calculated.
 */
template <size_t NumRegions, size_t ClosureOrder>
ScalarType get_central_mom_by_raw(
    std::array<int, static_cast<size_t>(osir::InfectionState::Count) * NumRegions> index,
    Eigen::Ref<const Eigen::VectorX<ScalarType>> y,
    const MomentArray<static_cast<size_t>(osir::InfectionState::Count), NumRegions, ClosureOrder>& moments,
    ScalarType E_XS_rS_XI_rI_XR_rR)
{
    int order = std::accumulate(index.begin(), index.end(), 0);
    if (order != ClosureOrder) {
        mio::log_error("The order of the moment that should be approximated is unequal to the given closure order. "
                       "Moment order is {} and closure order is {}.",
                       order, ClosureOrder);
    }
    auto rS           = index[0];
    auto rI           = index[1];
    auto rR           = index[2];
    ScalarType mean_S = y[0];
    ScalarType mean_I = y[1];
    ScalarType mean_R = y[2];
    double M_rS_rI_rR = 0.;
    for (int lS = 0; lS <= rS; ++lS) {
        for (int lI = 0; lI <= rI; ++lI) {
            for (int lR = 0; lR <= rR; ++lR) {
                double E_XS_lS_XI_lI_XR_lR = int(ClosureOrder) > (lS + lI + lR)
                                                 ? get_E_XS_lS_XI_lI_XR_lR(lS, lI, lR, y, moments)
                                                 : E_XS_rS_XI_rI_XR_rR;
                M_rS_rI_rR += boost::math::binomial_coefficient<double>(rS, lS) * std::pow(-1 * mean_S, rS - lS) *
                              boost::math::binomial_coefficient<double>(rI, lI) * std::pow(-1 * mean_I, rI - lI) *
                              boost::math::binomial_coefficient<double>(rR, lR) * std::pow(-1 * mean_R, rR - lR) *
                              E_XS_lS_XI_lI_XR_lR;
            }
        }
    }
    return M_rS_rI_rR;
}

/**
 * @brief Implements a closure for a given moment assuming that the underlying random variables are lognormally distributed. Closure only works for one region.
 * @tparam NumRegions Number of regions.
 * @tparam ClosureOrder Order to be closed.
 * @param[in] index Index of the central moment that should be approximated.
 * @param[in] y Current value of expected value and all moments.
 * @param[in] moments Moment array. Is only used to get the correct flat index of a moment in y.
 */
template <size_t NumRegions, size_t ClosureOrder>
ScalarType lognormal_closure(
    std::array<int, static_cast<size_t>(osir::InfectionState::Count) * NumRegions> index,
    Eigen::Ref<const Eigen::VectorX<ScalarType>> y,
    const MomentArray<static_cast<size_t>(osir::InfectionState::Count), NumRegions, ClosureOrder>& moments)
{
    int order = std::accumulate(index.begin(), index.end(), 0);
    if (order != ClosureOrder) {
        mio::log_error("The order of the moment that should be approximated is unequal to the given closure order. "
                       "Moment order is {} and closure order is {}.",
                       order, ClosureOrder);
    }
    if (NumRegions > 1) {
        mio::unused(index, y, moments);
        return 0;
    }
    auto rS                   = index[0];
    auto rI                   = index[1];
    auto rR                   = index[2];
    ScalarType mean_S         = y[0];
    ScalarType mean_I         = y[1];
    ScalarType mean_R         = y[2];
    ScalarType var_S          = y[moments.flatten_index({2, 0, 0}) + static_cast<size_t>(osir::InfectionState::Count)];
    ScalarType var_I          = y[moments.flatten_index({0, 2, 0}) + static_cast<size_t>(osir::InfectionState::Count)];
    ScalarType var_R          = y[moments.flatten_index({0, 0, 2}) + static_cast<size_t>(osir::InfectionState::Count)];
    ScalarType cov_SI         = y[moments.flatten_index({1, 1, 0}) + static_cast<size_t>(osir::InfectionState::Count)];
    ScalarType cov_SR         = y[moments.flatten_index({1, 0, 1}) + static_cast<size_t>(osir::InfectionState::Count)];
    ScalarType cov_IR         = y[moments.flatten_index({0, 1, 1}) + static_cast<size_t>(osir::InfectionState::Count)];
    ScalarType cov_SI_fac     = mean_S * mean_I != 0 ? std::pow(1 + cov_SI / (mean_S * mean_I), rS * rI) : 1.;
    ScalarType cov_SR_fac     = mean_S * mean_R != 0 ? std::pow(1 + cov_SR / (mean_S * mean_R), rS * rR) : 1.;
    ScalarType cov_IR_fac     = mean_I * mean_R != 0 ? std::pow(1 + cov_IR / (mean_I * mean_R), rI * rR) : 1.;
    ScalarType mean_fac_S     = (mean_S == 0.) && (2 * rS - rS * rS < 0) ? 0. : std::pow(mean_S, 2 * rS - rS * rS);
    ScalarType mean_fac_I     = (mean_I == 0.) && (2 * rI - rI * rI < 0) ? 0. : std::pow(mean_I, 2 * rI - rI * rI);
    ScalarType mean_fac_R     = (mean_R == 0.) && (2 * rR - rR * rR < 0) ? 0. : std::pow(mean_R, 2 * rR - rR * rR);
    ScalarType mean_var_fac_S = std::sqrt(var_S + mean_S * mean_S) == 0. && (rS * rS - rS < 0)
                                    ? 0.
                                    : std::pow(std::sqrt(var_S + mean_S * mean_S), rS * rS - rS);
    ScalarType mean_var_fac_I = std::sqrt(var_I + mean_I * mean_I) == 0. && (rI * rI - rI < 0)
                                    ? 0.
                                    : std::pow(std::sqrt(var_I + mean_I * mean_I), rI * rI - rI);
    ScalarType mean_var_fac_R = std::sqrt(var_R + mean_R * mean_R) == 0. && (rR * rR - rR < 0)
                                    ? 0.
                                    : std::pow(std::sqrt(var_R + mean_R * mean_R), rR * rR - rR);
    // Raw moment (r_S, r_I, r_R) which is approximated: E[X_S^r_S * X_I^r_I * X_R^r_R]
    ScalarType E_XS_rS_XI_rI_XR_rR = mean_fac_S * mean_var_fac_S * mean_fac_I * mean_var_fac_I * mean_fac_R *
                                     mean_var_fac_R * cov_SI_fac * cov_SR_fac * cov_IR_fac;
    // std::exp(rS * get_log_mu(var_S, mean_S) + rI * get_log_mu(var_I, mean_I) + rR * get_log_mu(var_R, mean_R) +
    //          0.5 * (rS * rS * std::pow(get_log_sigma(var_S, mean_S), 2.0) +
    //                 rI * rI * std::pow(get_log_sigma(var_I, mean_I), 2.0) +
    //                 rR * rR * std::pow(get_log_sigma(var_R, mean_R), 2.0)) +
    //          rS * rI * get_log_cov(cov_SI, mean_S, mean_I) + rS * rR * get_log_cov(cov_SR, mean_S, mean_R) +
    //          rI * rR * get_log_cov(cov_IR, mean_I, mean_R));

    return get_central_mom_by_raw(index, y, moments, E_XS_rS_XI_rI_XR_rR);
}

/**
 * @brief Implements a closure for a given moment assuming that the underlying random variables are lognormally distributed and assuming a mixture distribution with zero inflation for I.
 * @tparam NumRegions Number of regions.
 * @tparam ClosureOrder Order to be closed.
 * @param[in] index Index of the central moment that should be approximated.
 * @param[in] y Current value of expected value and all moments.
 * @param[in] moments Moment array. Is only used to get the correct flat index of a moment in y.
 */
template <size_t NumRegions, size_t ClosureOrder>
ScalarType lognormal_zero_inflation_closure(
    std::array<int, static_cast<size_t>(osir::InfectionState::Count) * NumRegions> index,
    Eigen::Ref<const Eigen::VectorX<ScalarType>> y,
    const MomentArray<static_cast<size_t>(osir::InfectionState::Count), NumRegions, ClosureOrder>& moments)
{
    int order = std::accumulate(index.begin(), index.end(), 0);
    if (order != ClosureOrder) {
        mio::log_error("The order of the moment that should be approximated is unequal to the given closure order. "
                       "Moment order is {} and closure order is {}.",
                       order, ClosureOrder);
    }
    if (NumRegions > 1) {
        mio::unused(index, y, moments);
        return 0;
    }
    auto rS               = index[0];
    auto rI               = index[1];
    auto rR               = index[2];
    ScalarType mean_S     = y[0];
    ScalarType mean_I     = y[1];
    ScalarType mean_R     = y[2];
    ScalarType var_S      = y[moments.flatten_index({2, 0, 0}) + static_cast<size_t>(osir::InfectionState::Count)];
    ScalarType var_I      = y[moments.flatten_index({0, 2, 0}) + static_cast<size_t>(osir::InfectionState::Count)];
    ScalarType var_R      = y[moments.flatten_index({0, 0, 2}) + static_cast<size_t>(osir::InfectionState::Count)];
    ScalarType cov_SI     = y[moments.flatten_index({1, 1, 0}) + static_cast<size_t>(osir::InfectionState::Count)];
    ScalarType cov_SR     = y[moments.flatten_index({1, 0, 1}) + static_cast<size_t>(osir::InfectionState::Count)];
    ScalarType cov_IR     = y[moments.flatten_index({0, 1, 1}) + static_cast<size_t>(osir::InfectionState::Count)];
    ScalarType skewness_S = y[moments.flatten_index({3, 0, 0}) + static_cast<size_t>(osir::InfectionState::Count)];
    ScalarType skewness_I = y[moments.flatten_index({0, 3, 0}) + static_cast<size_t>(osir::InfectionState::Count)];
    ScalarType skewness_R = y[moments.flatten_index({0, 0, 3}) + static_cast<size_t>(osir::InfectionState::Count)];
    ScalarType p          = 1. - (std::pow(mean_S, 6.) * (1. + 3. * var_S) + std::pow(mean_S, 3.) * skewness_S) /
                            std::pow(var_S + std::pow(mean_S, 2.), 3.);
    ScalarType cov_SI_fac = mean_S * mean_I != 0 ? std::pow(1 + cov_SI / (mean_S * mean_I), rS * rI) : 1.;
    ScalarType cov_SR_fac = mean_S * mean_R != 0 ? std::pow(1 + cov_SR / (mean_S * mean_R), rS * rR) : 1.;
    ScalarType cov_IR_fac = mean_I * mean_R != 0 ? std::pow(1 + cov_IR / (mean_I * mean_R), rI * rR) : 1.;
    ScalarType mean_fac_S = (mean_S == 0.) && (rS + rS * rS < 0) ? 0. : std::pow(std::sqrt(mean_S), rS + rS * rS);
    ScalarType mean_var_fac_S =
        (mean_S * mean_S + var_S == 0.) && (rS - rS * rS < 0) ? 0. : std::pow(mean_S * mean_S + var_S, rS - rS * rS);
    ScalarType mean_var_skew_fac_S =
        (mean_S * mean_S * mean_S + 3 * mean_S * var_S + skewness_S == 0.) && (rS - rS * rS < 0)
            ? 0.
            : std::pow(mean_S * mean_S * mean_S + 3 * mean_S * var_S + skewness_S, rS - rS * rS);
    ScalarType mean_fac_I = (mean_I == 0.) && (rI + rI * rI < 0) ? 0. : std::pow(std::sqrt(mean_I), rI + rI * rI);
    ScalarType mean_var_fac_I =
        (mean_I * mean_I + var_I == 0.) && (rI - rI * rI < 0) ? 0. : std::pow(mean_I * mean_I + var_I, rI - rI * rI);
    ScalarType mean_var_skew_fac_I =
        (mean_I * mean_I * mean_I + 3 * mean_I * var_I + skewness_I == 0.) && (rI - rI * rI < 0)
            ? 0.
            : std::pow(mean_I * mean_I * mean_I + 3 * mean_I * var_I + skewness_I, rI - rI * rI);
    ScalarType mean_fac_R = (mean_R == 0.) && (rR + rR * rR < 0) ? 0. : std::pow(std::sqrt(mean_R), rR + rR * rR);
    ScalarType mean_var_fac_R =
        (mean_R * mean_R + var_R == 0.) && (rR - rR * rR < 0) ? 0. : std::pow(mean_R * mean_R + var_R, rR - rR * rR);
    ScalarType mean_var_skew_fac_R =
        (mean_R * mean_R * mean_R + 3 * mean_R * var_R + skewness_R == 0.) && (rR - rR * rR < 0)
            ? 0.
            : std::pow(mean_R * mean_R * mean_R + 3 * mean_R * var_R + skewness_R, rR - rR * rR);

    // Raw moment (r_S, r_I, r_R) which is approximated: E[X_S^r_S * X_I^r_I * X_R^r_R]
    ScalarType E_XS_rS_XI_rI_XR_rR = std::pow(1 - p, 1 - rS - rI - rR + rS * rI + rS * rR + rI * rR) * mean_fac_S *
                                     mean_var_fac_S * mean_var_skew_fac_S * mean_fac_I * mean_var_fac_I *
                                     mean_var_skew_fac_I * mean_fac_R * mean_var_fac_R * mean_var_skew_fac_R *
                                     cov_SI_fac * cov_SR_fac * cov_IR_fac;
    for (auto& i : index) {
        std::cout << " " << i;
    }
    std::cout << ": " << E_XS_rS_XI_rI_XR_rR << std::endl;

    return get_central_mom_by_raw(index, y, moments, E_XS_rS_XI_rI_XR_rR);
}

/**
 * @brief Implements a closure for a given moment assuming independence of the underlying random variables. Closure only works for one region.
 * @tparam NumRegions Number of regions.
 * @tparam ClosureOrder Order to be closed.
 * @param[in] index Index of the central moment that should be approximated.
 * @param[in] y Current value of expected value and all moments.
 * @param[in] moments Moment array. Is only used to get the correct flat index of a moment in y.
 */
template <size_t NumRegions, size_t ClosureOrder>
ScalarType pairapprox_closure(
    std::array<int, static_cast<size_t>(osir::InfectionState::Count) * NumRegions> index,
    Eigen::Ref<const Eigen::VectorX<ScalarType>> y,
    const MomentArray<static_cast<size_t>(osir::InfectionState::Count), NumRegions, ClosureOrder>& moments)
{
    int order = std::accumulate(index.begin(), index.end(), 0);
    if (order != ClosureOrder) {
        mio::log_error("The order of the moment that should be approximated is unequal to the given closure order. "
                       "Moment order is {} and closure order is {}.",
                       order, ClosureOrder);
    }
    if (NumRegions > 1) {
        mio::unused(index, y, moments);
        return 0;
    }
    auto rS           = index[0];
    auto rI           = index[1];
    auto rR           = index[2];
    ScalarType mean_S = y[0];
    ScalarType mean_I = y[1];
    ScalarType mean_R = y[2];
    // Raw moment (r_S, r_I, r_R) which is approximated: E[X_S^r_S * X_I^r_I * X_R^r_R]
    ScalarType E_XS_rS_XI_rI_XR_rR = std::pow(mean_S, rS) * std::pow(mean_I, rI) * std::pow(mean_R, rR);

    return get_central_mom_by_raw(index, y, moments, E_XS_rS_XI_rI_XR_rR);
}

} // namespace smm_moments
} // namespace mio

#endif // MOMENTS_CLOSURE_FUNCTIONS_H
