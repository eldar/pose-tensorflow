#pragma once
#ifndef ANDRES_GRAPH_PATHS_HXX
#define ANDRES_GRAPH_PATHS_HXX

#include <cstddef>
#include <utility> // std::pair

#include "andres/graph/graph.hxx" // DefaultSubgraphMask

namespace andres {
namespace graph {

/// Search a path for a chord.
///
/// \param graph Graph.
/// \param begin Iterator to the beginning of the sequence of nodes on the path.
/// \param end Iterator to the end of the sequence of nodes on the path.
/// \param ignoreEdgeBetweenFirstAndLast Flag.
///
template<class GRAPH, class ITERATOR>
inline std::pair<bool, std::size_t>
findChord(
    const GRAPH& graph,
    ITERATOR begin,
    ITERATOR end,
    const bool ignoreEdgeBetweenFirstAndLast = false
) {
    return findChord(graph, DefaultSubgraphMask<>(), begin, end, 
        ignoreEdgeBetweenFirstAndLast);
}

/// Search a path for a chord.
///
/// \param graph Graph.
/// \param mask A subgraph mask such as DefaultSubgraphMask.
/// \param begin Iterator to the beginning of the sequence of nodes on the path.
/// \param end Iterator to the end of the sequence of nodes on the path.
/// \param ignoreEdgeBetweenFirstAndLast Flag.
///
template<class GRAPH, class SUBGRAPH_MASK, class ITERATOR>
inline std::pair<bool, std::size_t>
findChord(
    const GRAPH& graph,
    const SUBGRAPH_MASK& mask,
    ITERATOR begin,
    ITERATOR end,
    const bool ignoreEdgeBetweenFirstAndLast = false
) {
    for(ITERATOR it = begin; it != end - 1; ++it) 
    for(ITERATOR it2 = it + 2; it2 != end; ++it2) {
        if(ignoreEdgeBetweenFirstAndLast && it == begin && it2 == end - 1) {
            continue;
        }
        std::pair<bool, std::size_t> p = graph.findEdge(*it, *it2);
        if(p.first && mask.edge(p.second)) {
            return p;
        }
    }
    return std::pair<bool, std::size_t>(false, 0);
}

} // namespace graph
} // namespace andres

#endif // #ifndef ANDRES_GRAPH_PATHS_HXX
