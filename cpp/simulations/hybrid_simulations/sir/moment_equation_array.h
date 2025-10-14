#include "memilio/config.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/utils/index.h"
#include "memilio/utils/time_series.h"
#include "ode_sir/infection_state.h"
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Meta.h>
#include <cassert>
#include <cstddef>
#include <vector>

/** @brief Structure to define a custom index type for the moment equation array; is just a typesafe wrapper for size_t.
 *  This allows us to create a custom index type that can be used with CustomIndexArray.
 */
struct MyIndex : public mio::Index<MyIndex> {
    MyIndex(int val)
        : mio::Index<MyIndex>(val)
    {
    }
};

/**
 * @brief A struct used to create a custom index array for moment equations.
 * The moments have one index per infection state i.e. one dimension per infection state is needed, and the order of the moment is defined by the template parameter Order.
 * @tparam Order The maximum order of the moment equations.
 * @tparam Is A sequence of length InfectionState::Count.
 */
template <int Order, std::size_t... Is>
struct repeat_index {
    using type = mio::CustomIndexArray<ScalarType, decltype((void)Is, MyIndex{Order})...>;
    // (void)Is just to silence unused warning — creates N Index types
};

// helper to expand the index sequence
template <std::size_t N, int Order, std::size_t... Is>
auto make_repeat_index(std::index_sequence<Is...>) -> typename repeat_index<Order, Is...>::type;

/**
 * @brief A helper struct to create a custom index array for moment equations.
 * It uses the repeat_index struct to create a sequence of indices for each infection state.
 * @tparam N The number of infection states, InfectionState::Count.
 * @tparam Order The maximum order of the moment equations.
 */
template <std::size_t N, int Order>
struct MakeArray {
    using type = decltype(make_repeat_index<N, Order>(std::make_index_sequence<N>{}));
};

/**
 * @brief A class template for an array of moment equations and equation for the expected value.
 * @tparam InfectionState The infection state enum.
 * @tparam Order The maximum order of the moment equations.
 */
template <class InfectionState, size_t NumRegions, int Order>
class MomentEquationArray
{
    // Number of infection states and regions defining the size of the expected values array, the number of moment indices and hence the number of dimensions in the CustomIndexArray.
    static constexpr std::size_t count = static_cast<std::size_t>(InfectionState::Count) * NumRegions;

    // Build an MType::Index with N copies of MyIndex(Order)
    template <class Index, std::size_t... Is>
    constexpr Index make_index_impl(int order, std::index_sequence<Is...>)
    {
        return Index{((void)Is, MyIndex(order))...};
    }

    template <std::size_t N, class Index>
    constexpr Index make_index(int order)
    {
        return make_index_impl<Index>(order, std::make_index_sequence<N>{});
    }

public:
    // CustomIndexArray for moment equations, with one index per infection state and the order of the moment defined by the template parameter Order.
    using MType = typename MakeArray<count, Order>::type;

    MomentEquationArray(std::vector<double> initial_populations)
        : m_moments(make_index<count, typename MType::Index>(Order), 0.0)
    {
        assert(Order >= 2);
        assert(initial_populations.size() == count);
        m_moments[make_index<count, typename MType::Index>(0)] = 1.0; // Set the zeroth moment to 1
        m_means = Eigen::Map<const Eigen::Array<ScalarType, Eigen::Dynamic, 1>>(initial_populations.data(),
                                                                                initial_populations.size());
    }

    /**
     * @brief Return a vector with the expected values and the moments.
     */
    Eigen::VectorXd get_values() const
    {
        Eigen::VectorXd combined(m_means.size() + m_moments.array().size());
        combined << m_means, m_moments.array();
        return combined;
    }

    /**
     * @brief Get the flat index for a given moment equation index.
     */
    size_t get_moment_index(MType::Index index) const
    {
        size_t flat_index = m_moments.get_flat_index(index);
        // Add the offset for the expected values, which are stored before the moment equations.
        return static_cast<size_t>(InfectionState::Count) * NumRegions + flat_index;
    }

private:
    MType m_moments; // CustomIndexArray for moment equations
    Eigen::Array<ScalarType, Eigen::Dynamic, 1> m_means; // Array for expected values equations
};

#include <boost/math/special_functions/binomial.hpp>

template double boost::math::binomial_coefficient<double>(unsigned, unsigned);

namespace SIR
{

double get_moment(const Eigen::Matrix<ScalarType, Eigen::Dynamic, 3>& values, const std::vector<int>& indices)
{
    // values has 3 columns: S, I, R
    Eigen::Matrix<ScalarType, 1, 3> means = values.colwise().mean();
    double moment                         = 0.0;
    for (int i = 0; i < values.rows(); ++i) {
        moment += std::pow(values(i, 0) - means(0), indices[0]) * std::pow(values(i, 1) - means(1), indices[1]) *
                  std::pow(values(i, 2) - means(2), indices[2]);
    }
    moment /= static_cast<double>(values.rows());
    return moment;
}

template <int Order>
class Moments
{
public:
    Moments(std::vector<double> initial_populations)
        : values(initial_populations)
    {
    }

    void get_derivatives(Eigen::Ref<const Eigen::VectorX<ScalarType>> y, ScalarType /*t*/,
                         Eigen::Ref<Eigen::VectorX<ScalarType>> dydt, int closure_order, double lambda,
                         double gamma) const
    {
        size_t mu_S = static_cast<size_t>(mio::osir::InfectionState::Susceptible);
        size_t mu_I = static_cast<size_t>(mio::osir::InfectionState::Infected);
        size_t mu_R = static_cast<size_t>(mio::osir::InfectionState::Recovered);
        assert(closure_order >= 2);
        double M_110 = closure_order > 2 ? y[values.get_moment_index({MyIndex(1), MyIndex(1), MyIndex(0)})] : 0;
        dydt[mu_S]   = -lambda * y[mu_S] * y[mu_I] - lambda * M_110;
        dydt[mu_I]   = lambda * y[mu_S] * y[mu_I] - gamma * y[mu_I] + lambda * M_110;
        dydt[mu_R]   = gamma * y[mu_I];
        if (closure_order > 2) {
            for (int i_S = 0; i_S < closure_order; i_S++) {
                for (int i_I = 0; i_I < closure_order; i_I++) {
                    for (int i_R = 0; i_R < closure_order; i_R++) {
                        if (i_S + i_I + i_R == 0 || i_S + i_I + i_R == 1 || i_S + i_I + i_R >= closure_order) {
                            continue;
                        }
                        else {
                            size_t index = values.get_moment_index({MyIndex(i_S), MyIndex(i_I), MyIndex(i_R)});
                            dydt[index]  = 0;
                            for (int l_S = 0; l_S <= i_S; l_S++) {
                                for (int l_I = 0; l_I <= i_I; l_I++) {
                                    if (!(l_S + l_I == i_S + i_I)) {
                                        if (!(l_S + l_I + i_R == 1)) {
                                            double M_lS_lI_iR = (l_S + l_I + i_R) >= closure_order
                                                                    ? 0.
                                                                    : y[values.get_moment_index(
                                                                          {MyIndex(l_S), MyIndex(l_I), MyIndex(i_R)})];
                                            if (l_S + l_I + i_R == 0) {
                                                M_lS_lI_iR = 1.;
                                            }
                                            dydt[index] += lambda * y[mu_S] * y[mu_I] *
                                                           boost::math::binomial_coefficient<double>(i_S, l_S) *
                                                           boost::math::binomial_coefficient<double>(i_I, l_I) *
                                                           std::pow(-1, i_S - l_S) * M_lS_lI_iR;
                                        }
                                        if (!(l_S + 1 + l_I + i_R == 1)) {
                                            double M_lSp1_lI_iR =
                                                (l_S + 1 + l_I + i_R) >= closure_order
                                                    ? 0.
                                                    : y[values.get_moment_index(
                                                          {MyIndex(l_S + 1), MyIndex(l_I), MyIndex(i_R)})];
                                            if (l_S + 1 + l_I + i_R == 0) {
                                                M_lSp1_lI_iR = 1.;
                                            }
                                            dydt[index] += lambda * y[mu_I] *
                                                           boost::math::binomial_coefficient<double>(i_S, l_S) *
                                                           boost::math::binomial_coefficient<double>(i_I, l_I) *
                                                           std::pow(-1, i_S - l_S) * M_lSp1_lI_iR;
                                        }
                                        if (!(l_S + l_I + 1 + i_R == 1)) {
                                            double M_lS_lIp1_iR =
                                                (l_S + l_I + 1 + i_R) >= closure_order
                                                    ? 0.
                                                    : y[values.get_moment_index(
                                                          {MyIndex(l_S), MyIndex(l_I + 1), MyIndex(i_R)})];
                                            if (l_S + l_I + 1 + i_R == 0) {
                                                M_lS_lIp1_iR = 1.;
                                            }
                                            dydt[index] += lambda * y[mu_S] *
                                                           boost::math::binomial_coefficient<double>(i_S, l_S) *
                                                           boost::math::binomial_coefficient<double>(i_I, l_I) *
                                                           std::pow(-1, i_S - l_S) * M_lS_lIp1_iR;
                                        }
                                        if (!(l_S + 1 + l_I + 1 + i_R == 1)) {
                                            double M_lSp1_lIp1_iR =
                                                (l_S + 1 + l_I + 1 + i_R) >= closure_order
                                                    ? 0.
                                                    : y[values.get_moment_index(
                                                          {MyIndex(l_S + 1), MyIndex(l_I + 1), MyIndex(i_R)})];
                                            if (l_S + 1 + l_I + 1 + i_R == 0) {
                                                M_lSp1_lIp1_iR = 1.;
                                            }
                                            dydt[index] += lambda *
                                                           boost::math::binomial_coefficient<double>(i_S, l_S) *
                                                           boost::math::binomial_coefficient<double>(i_I, l_I) *
                                                           std::pow(-1, i_S - l_S) * M_lSp1_lIp1_iR;
                                        }
                                    }
                                }
                            }
                            for (int l_I = 0; l_I <= i_I; l_I++) {
                                for (int l_R = 0; l_R <= i_R; l_R++) {
                                    if (!(l_I + l_R == i_I + i_R)) {
                                        if (!(i_S + l_I + l_R == 1)) {
                                            double M_iS_lI_lR = (i_S + l_I + l_R) >= closure_order
                                                                    ? 0.
                                                                    : y[values.get_moment_index(
                                                                          {MyIndex(i_S), MyIndex(l_I), MyIndex(l_R)})];
                                            if (i_S + l_I + l_R == 0) {
                                                M_iS_lI_lR = 1.;
                                            }
                                            dydt[index] += gamma * y[mu_I] *
                                                           boost::math::binomial_coefficient<double>(i_I, l_I) *
                                                           boost::math::binomial_coefficient<double>(i_R, l_R) *
                                                           std::pow(-1, i_I - l_I) * M_iS_lI_lR;
                                        }

                                        if (!(i_S + l_I + 1 + l_R == 1)) {
                                            double M_iS_lIp1_lR =
                                                (i_S + l_I + 1 + l_R) >= closure_order
                                                    ? 0.
                                                    : y[values.get_moment_index(
                                                          {MyIndex(i_S), MyIndex(l_I + 1), MyIndex(l_R)})];
                                            if (i_S + l_I + 1 + l_R == 0) {
                                                M_iS_lIp1_lR = 1.;
                                            }
                                            dydt[index] += gamma * boost::math::binomial_coefficient<double>(i_I, l_I) *
                                                           boost::math::binomial_coefficient<double>(i_R, l_R) *
                                                           std::pow(-1, i_I - l_I) * M_iS_lIp1_lR;
                                        }
                                    }
                                }
                            }
                            if (!(i_S - 1 + i_I + i_R == 1)) {
                                double M_iSm1_iI_iR = 0.;
                                if (i_S > 0) {
                                    M_iSm1_iI_iR =
                                        y[values.get_moment_index({MyIndex(i_S - 1), MyIndex(i_I), MyIndex(i_R)})];
                                }
                                double M_iS_iIm1_iR = 0.;
                                if (i_I > 0) {
                                    M_iS_iIm1_iR =
                                        y[values.get_moment_index({MyIndex(i_S), MyIndex(i_I - 1), MyIndex(i_R)})];
                                }
                                double M_iS_iI_iRm1 = 0.;
                                if (i_R > 0) {
                                    M_iS_iI_iRm1 =
                                        y[values.get_moment_index({MyIndex(i_S), MyIndex(i_I), MyIndex(i_R - 1)})];
                                }
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

    mio::TimeSeries<ScalarType> calculate_moments_from_sim(const std::vector<mio::TimeSeries<ScalarType>>& sim_results)
    {
        mio::TimeSeries<double> moments_ts(values.get_values().size());
        for (int t = 0; t < sim_results[0].get_num_time_points(); ++t) {
            Eigen::VectorXd moment_values = Eigen::VectorXd::Zero(values.get_values().size());
            Eigen::Matrix<ScalarType, Eigen::Dynamic, 3> result =
                Eigen::Matrix<ScalarType, Eigen::Dynamic, 3>::Zero(sim_results.size(), 3);
            for (size_t i = 0; i < sim_results.size(); i++) {
                result(i, 0) = sim_results[i].get_value(t)[static_cast<size_t>(mio::osir::InfectionState::Susceptible)];
                result(i, 1) = sim_results[i].get_value(t)[static_cast<size_t>(mio::osir::InfectionState::Infected)];
                result(i, 2) = sim_results[i].get_value(t)[static_cast<size_t>(mio::osir::InfectionState::Recovered)];
            }
            moment_values[0] = result.col(0).mean(); // muS
            moment_values[1] = result.col(1).mean(); // muI
            moment_values[2] = result.col(2).mean(); // muR
            for (int i_S = 0; i_S < Order; ++i_S) {
                for (int i_I = 0; i_I < Order; ++i_I) {
                    for (int i_R = 0; i_R < Order; ++i_R) {
                        if (i_S + i_I + i_R >= Order) {
                            continue;
                        }
                        std::vector<int> indices{i_S, i_I, i_R};
                        size_t index         = values.get_moment_index({MyIndex(i_S), MyIndex(i_I), MyIndex(i_R)});
                        moment_values[index] = get_moment(result, indices);
                    }
                }
            }
            moments_ts.add_time_point(sim_results[0].get_time(t), std::move(moment_values));
        }
        return moments_ts;
    }

    std::pair<mio::TimeSeries<double>, std::vector<std::string>>
    get_moment_ts_and_names(const mio::TimeSeries<double>& sim_result_ts, size_t closing_order)
    {
        int num_moments = 3; // muS, muI, muR
        for (int m = 0; m < (int)closing_order; ++m) {
            num_moments += boost::math::binomial_coefficient<double>(m + 2, 2);
        }
        // Map having the index of a moment it has in sim_result_ts as key and the name of the moment as value
        std::map<size_t, std::string> index_to_names;
        index_to_names[0] = "muS";
        index_to_names[1] = "muI";
        index_to_names[2] = "muR";
        // Timeseries corresponding only expected values and moments up to closing order-1
        mio::TimeSeries<double> moment_ts =
            mio::TimeSeries<double>::zero(sim_result_ts.get_num_time_points(), num_moments);
        // Get moments names and corresponding indices in sim_result_ts
        for (int iS = 0; iS < (int)closing_order; ++iS) {
            for (int iI = 0; iI < (int)closing_order; ++iI) {
                for (int iR = 0; iR < (int)closing_order; ++iR) {
                    if (iS + iI + iR >= int(closing_order)) {
                        continue;
                    }
                    size_t index          = values.get_moment_index({MyIndex(iS), MyIndex(iI), MyIndex(iR)});
                    index_to_names[index] = "M" + std::to_string(iS) + std::to_string(iI) + std::to_string(iR);
                }
            }
        }

        // Fill the new timeseries with only the expected values and moments needed
        for (auto t = Eigen::Index(0); t < sim_result_ts.get_num_time_points(); ++t) {
            size_t index          = 0;
            moment_ts.get_time(t) = sim_result_ts.get_time(t);
            for (const auto& [key, value] : index_to_names) {
                moment_ts.get_value(t)[index] = sim_result_ts.get_value(t)[key];
                ++index;
            }
        }

        // Get the names of the moments in the correct order
        std::vector<std::string> names;
        names.reserve(index_to_names.size());

        for (const auto& [key, value] : index_to_names) {
            names.push_back(value);
        }
        return std::make_pair(moment_ts, names);
    }

    MomentEquationArray<mio::osir::InfectionState, 1, Order> values;
};

} // namespace SIR
