#include "memilio/config.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/utils/index.h"
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
template <class InfectionState, int Order>
class MomentEquationArray
{
    // Number of infection states defining the size of the expected values array, the number of moment indices and hence the number of dimensions in the CustomIndexArray.
    static constexpr std::size_t count = static_cast<std::size_t>(InfectionState::Count);
    // CustomIndexArray for moment equations, with one index per infection state and the order of the moment defined by the template parameter Order.
    using MType = typename MakeArray<count, Order>::type;

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
    MomentEquationArray(std::vector<double> initial_populations)
        : m_moments(make_index<count, typename MType::Index>(Order), 0.0)
    {
        assert(initial_populations.size() == count);
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
        return int(InfectionState::Count) + m_moments.get_flat_index(index);
    }

private:
    MType m_moments; // CustomIndexArray for moment equations
    Eigen::Array<ScalarType, Eigen::Dynamic, 1> m_means; // Array for expected values equations
};
