#pragma once
#ifndef ANDRES_GRAPH_EDGE_VALUE_HXX
#define ANDRES_GRAPH_EDGE_VALUE_HXX

namespace andres {
namespace graph {

/// Return 1 for every edge.
template<class T>
struct UnitEdgeValueIterator {
    typedef std::ptrdiff_t difference_type;
    typedef T value_type;
    typedef value_type* pointer;
    typedef value_type& reference;
    typedef std::random_access_iterator_tag iterator_category;

    value_type operator[](const std::size_t) const
        { return 1; }
};

} // namespace graph
} // namespace andres

#endif // #ifndef ANDRES_GRAPH_EDGE_VALUE_HXX
