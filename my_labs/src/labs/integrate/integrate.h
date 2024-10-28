#ifndef INTEGRATE_H
#define INTEGRATE_H

#include <array>
#include <type_traits>
#include <cstdint>
#include <cmath>
using u32 = std::uint32_t;

template<typename A>
struct ArgumentGetter;

template<typename R, typename Arg>
struct ArgumentGetter<R(Arg)> {
    using Argument = Arg;
};

template<typename T>
using Dif = decltype(std::declval<T>() - std::declval<T>());

template <typename Callable>
struct Node
{
    typename ArgumentGetter<Callable>::Argument point;
    double weight;
};

template <typename Callable, u32 N>
struct Nodes
{
    static std::array<Node<Callable>, N> nodes;
};

template <typename Callable>
struct Nodes<Callable, 1>
{
    static constexpr std::array<Node<Callable>, 1> nodes = {Node<Callable>{0.5, 1.}};
};

template <typename Callable>
struct Nodes<Callable, 2>
{
    static constexpr std::array<Node<Callable>, 2> nodes =
    {
        Node<Callable>{0.21132486540518711774, 0.5},
        Node<Callable>{0.78867513459481288225, 0.5}
    };
};

template <typename Callable>
struct Nodes<Callable, 3>
{
    static constexpr std::array<Node<Callable>, 3> nodes =
    {
        Node<Callable>{0.11270166537925831148, 0.27777777777777777777},
        Node<Callable>{0.50000000000000000000, 0.44444444444444444444},
        Node<Callable>{0.88729833462074168851, 0.27777777777777777777},
    };
};

template <typename Callable>
struct Nodes<Callable, 4>
{
    static constexpr std::array<Node<Callable>, 4> nodes =
    {
        Node<Callable>{0.06943184420297371238, 0.17392742256872692868},
        Node<Callable>{0.33000947820757186759, 0.32607257743127307131},
        Node<Callable>{0.66999052179242813240, 0.32607257743127307131},
        Node<Callable>{0.93056815579702628761, 0.17392742256872692868}
    };
};

template <typename Callable>
struct Nodes<Callable, 5>
{
    static constexpr std::array<Node<Callable>, 5> nodes =
    {
        Node<Callable>{0.04691007703066800360, 0.11846344252809454375},
        Node<Callable>{0.23076534494715845448, 0.23931433524968323402},
        Node<Callable>{0.50000000000000000000, 0.28444444444444444444},
        Node<Callable>{0.76923465505284154551, 0.23931433524968323402},
        Node<Callable>{0.95308992296933199639, 0.11846344252809454375}
    };
};

template <typename Callable>
struct Nodes<Callable, 6>
{
    static constexpr std::array<Node<Callable>, 6> nodes =
    {
        Node<Callable>{0.03376524289842398609, 0.08566224618958517252},
        Node<Callable>{0.16939530676686774316, 0.18038078652406930378},
        Node<Callable>{0.38069040695840154568, 0.23395696728634552369},
        Node<Callable>{0.61930959304159845431, 0.23395696728634552369},
        Node<Callable>{0.83060469323313225683, 0.18038078652406930378},
        Node<Callable>{0.96623475710157601390, 0.08566224618958517252}
    };
};

template <typename Callable>
struct Nodes<Callable, 7>
{
    static constexpr std::array<Node<Callable>, 7> nodes =
    {
        Node<Callable>{0.02544604382862073773, 0.06474248308443484663},
        Node<Callable>{0.12923440720030278006, 0.13985269574463833395},
        Node<Callable>{0.29707742431130141654, 0.19091502525255947247},
        Node<Callable>{0.50000000000000000000, 0.20897959183673469387},
        Node<Callable>{0.70292257568869858345, 0.19091502525255947247},
        Node<Callable>{0.87076559279969721993, 0.13985269574463833395},
        Node<Callable>{0.97455395617137926226, 0.06474248308443484663}
    };
};

template <typename Callable>
struct Nodes<Callable, 8>
{
    static constexpr std::array<Node<Callable>, 8> nodes =
    {
        Node<Callable>{0.01985507175123188415, 0.05061426814518812957},
        Node<Callable>{0.10166676129318663020, 0.11119051722668723527},
        Node<Callable>{0.23723379504183550709, 0.15685332293894364366},
        Node<Callable>{0.40828267875217509753, 0.18134189168918099148},
        Node<Callable>{0.59171732124782490246, 0.18134189168918099148},
        Node<Callable>{0.76276620495816449290, 0.15685332293894364366},
        Node<Callable>{0.89833323870681336979, 0.11119051722668723527},
        Node<Callable>{0.98014492824876811584, 0.05061426814518812957}
    };
};

template<u32 N, typename Callable>
decltype(auto) integrate(
    const Callable& func,
    const typename ArgumentGetter<Callable>::Argument& start,
    const typename ArgumentGetter<Callable>::Argument& end
    )
{
    Dif<typename ArgumentGetter<Callable>::Argument> h = end - start;
    double sum = 0;
    for (u32 i = 0; i < N; ++i)
    {
        sum += func(start + h * Nodes<Callable, N>::nodes[i].point) * Nodes<Callable, N>::nodes[i].weight;
    }
    return sum*h;
}

template<u32 N, typename Callable>
decltype(auto) integrate(
    const Callable& func,
    const typename ArgumentGetter<Callable>::Argument& start,
    const typename ArgumentGetter<Callable>::Argument& end,
    const Dif<typename ArgumentGetter<Callable>::Argument>& dx
    )
{
    double integral = 0;
    for (auto x = start; x < end; x += static_cast<typename ArgumentGetter<Callable>::Argument>(dx))
    {
        if (x + dx > end)
        {
            integral += integrate<N>(func, x, end);
            break;
        }
        integral += integrate<N>(func, x, x + dx);
    }
    return integral;
}

template<u32 N, typename Callable>
decltype(auto) integrateRichardsonExtrapolation(
    const Callable& func,
    typename ArgumentGetter<Callable>::Argument const &start,
    typename ArgumentGetter<Callable>::Argument const &end,
    const Dif<typename ArgumentGetter<Callable>::Argument>& dx
    )
{
    auto I     = integrate<N>(func, start, end, dx);
    auto Ihalf = integrate<N>(func, start, end, dx/2);

    double denom = std::pow(2, 2*N+1) - 1;
    auto delta = (I - Ihalf)/denom;
    return Ihalf - delta;
}

template<u32 N, typename Callable>
decltype(auto) integrate(
    const Callable& func,
    typename ArgumentGetter<Callable>::Argument const &start,
    typename ArgumentGetter<Callable>::Argument const &end,
    double const &err,
    u32 const &maxIt
    )
{
    Dif<typename ArgumentGetter<Callable>::Argument> h = end - start;
    auto I     = integrate<N>(func, start, end, h);
    auto Ihalf = integrate<N>(func, start, end, h/2);

    double denom = std::pow(2, 2*N+1) - 1;
    auto delta = (I - Ihalf)/denom;
    for (u32 i = 0; i < maxIt, std::abs(delta) > err; ++i)
    {
        h = h/2;
        I     = Ihalf;
        Ihalf = integrate<N>(func, start, end, h/2);
        delta = (I - Ihalf)/denom;
    }
    return Ihalf - delta;
}

#endif //INTEGRATE_H
