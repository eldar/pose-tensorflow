#pragma once
#ifndef ANDRES_GRAPH_SUBGRAPH_HXX
#define ANDRES_GRAPH_SUBGRAPH_HXX

namespace andres {
namespace graph {

/// An entire graph.
template<class T = std::size_t>
struct DefaultSubgraphMask {
    typedef T Value;

    bool vertex(const Value v) const
        { return true; }
    bool edge(const Value e) const
        { return true; }
};

} // namespace graph
} // namespace andres

#endif // #ifndef ANDRES_GRAPH_SUBGRAPH_HXX
