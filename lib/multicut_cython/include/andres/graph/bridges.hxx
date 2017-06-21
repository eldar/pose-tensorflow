#pragma once
#ifndef ANDRES_GRAPH_BRIDGES_HXX
#define ANDRES_GRAPH_BRIDGES_HXX

#include <stack>
#include <vector>
#include "subgraph.hxx"

namespace andres {
namespace graph {

template<class GRAPH>
struct BridgesBuffers
{
    typedef GRAPH GraphType;

    BridgesBuffers(const GraphType&);

    std::vector<std::size_t> depth_;
    std::vector<std::size_t> min_successor_depth_;
    std::vector<typename GraphType::VertexIterator> next_out_arc_;
    std::vector<int> parent_;
    std::vector<char> visited_;
};

template<typename GRAPH, typename ARRAY>
inline void
findBridges(
    const GRAPH&,
    ARRAY&
);

template<typename GRAPH, typename ARRAY>
inline void
findBridges(
    const GRAPH&,
    ARRAY&,
    BridgesBuffers<GRAPH>&
);

template<typename GRAPH, typename SUBGRAPH, typename ARRAY>
inline void
findBridges(
    const GRAPH&,
    const SUBGRAPH&,
    ARRAY&
);

template<typename GRAPH, typename SUBGRAPH, typename ARRAY>
inline void
findBridges(
    const GRAPH&,
    const SUBGRAPH&,
    ARRAY&,
    BridgesBuffers<GRAPH>&
);

template<typename GRAPH, typename SUBGRAPH, typename ARRAY>
inline void
findBridges(
    const GRAPH&,
    const SUBGRAPH&,
    std::size_t,
    ARRAY&
);

template<typename GRAPH, typename SUBGRAPH, typename ARRAY>
inline void
findBridges(
    const GRAPH&,
    const SUBGRAPH&,
    std::size_t,
    ARRAY&,
    BridgesBuffers<GRAPH>&
);

template<typename GRAPH>
BridgesBuffers<GRAPH>::BridgesBuffers(const GraphType& graph) :
    depth_(graph.numberOfVertices()),
    min_successor_depth_(graph.numberOfVertices()),
    next_out_arc_(graph.numberOfVertices()),
    parent_(graph.numberOfVertices()),
    visited_(graph.numberOfVertices())
{}

/// \brief Find, by Tarjan's algorithm, bridges (cut edges) in an undirected graph.
///
/// Tarjan, R. (1974). A note on finding the bridges of a graph.
/// Information Processing Letters 2(6):160-161.
///
/// Runtime complexity O(|E| + |V|).
/// 
/// \param graph An undirected graph.
/// \param is_bridge Array storing 1 for each edge-index if it is a bridge, or 0 otherwise.
///
template<typename GRAPH, typename ARRAY>
inline void
findBridges(
    const GRAPH& graph,
    ARRAY& is_bridge
)
{
    auto buffer = BridgesBuffers<GRAPH>(graph);
    findBridges(graph, DefaultSubgraphMask<>(), is_bridge, buffer);
}

/// \brief Find, by Tarjan's algorithm, bridges (cut edges) in an undirected graph.
///
/// Tarjan, R. (1974). A note on finding the bridges of a graph.
/// Information Processing Letters 2(6):160-161.
///
/// Runtime complexity O(|E| + |V|).
/// 
/// \param graph An undirected graph.
/// \param is_bridge Array storing 1 for each edge-index if it is a bridge, or 0 otherwise.
/// \param buffer A pre-allocated buffer object.
///
template<typename GRAPH, typename ARRAY>
inline void
findBridges(
    const GRAPH& graph,
    ARRAY& is_bridge,
    BridgesBuffers<GRAPH>& buffer
)
{
    findBridges(graph, DefaultSubgraphMask<>(), is_bridge, buffer);
}

/// \brief Find, by Tarjan's algorithm, bridges (cut edges) in a subgraph of an undirected graph.
///
/// Tarjan, R. (1974). A note on finding the bridges of a graph.
/// Information Processing Letters 2(6):160-161.
///
/// Runtime complexity O(|E| + |V|).
/// 
/// \param graph An undirected graph.
/// \param subgraph_mask Mask defining the subgraph.
/// \param is_bridge Array storing 1 for each edge-index if it is a bridge, or 0 otherwise.
///
template<typename GRAPH, typename SUBGRAPH, typename ARRAY>
inline void
findBridges(
    const GRAPH& graph,
    const SUBGRAPH& subgraph_mask,
    ARRAY& is_bridge
)
{
    auto buffer = BridgesBuffers<GRAPH>(graph);
    findBridges(graph, subgraph_mask, is_bridge, buffer);
}

/// \brief Find, by Tarjan's algorithm, bridges (cut edges) in a subgraph of an undirected graph.
///
/// Tarjan, R. (1974). A note on finding the bridges of a graph.
/// Information Processing Letters 2(6):160-161.
///
/// Runtime complexity O(|E| + |V|).
/// 
/// \param graph An undirected graph.
/// \param subgraph_mask Mask defining the subgraph.
/// \param is_bridge Array storing 1 for each edge-index if it is a bridge, or 0 otherwise.
/// \param buffer A pre-allocated buffer object.
///
template<typename GRAPH, typename SUBGRAPH, typename ARRAY>
inline void
findBridges(
    const GRAPH& graph,
    const SUBGRAPH& subgraph_mask,
    ARRAY& is_bridge,
    BridgesBuffers<GRAPH>& buffer
)
{
   std::fill(buffer.parent_.begin(), buffer.parent_.end(), -2);

    for (std::size_t i = 0; i < graph.numberOfVertices(); ++i)
        if (buffer.parent_[i] == -2 && subgraph_mask.vertex(i))
            findBridges(graph, subgraph_mask, i, is_bridge, buffer); 
}

/// \brief Find, by Tarjan's algorithm, bridges (cut edges) in a subgraph of an undirected graph containing a vertex.
///
/// Tarjan, R. (1974). A note on finding the bridges of a graph.
/// Information Processing Letters 2(6):160-161.
///
/// Runtime complexity O(|E| + |V|).
/// 
/// \param graph An undirected graph.
/// \param subgraph_mask Mask defining the subgraph.
/// \param starting_vertex A vertex inside the subgraph.
/// \param is_bridge Array storing 1 for each edge-index if it is a bridge, or 0 otherwise.
///
template<typename GRAPH, typename SUBGRAPH, typename ARRAY>
inline void
findBridges(
    const GRAPH& graph,
    const SUBGRAPH& subgraph_mask,
    std::size_t starting_vertex,
    ARRAY& is_bridge
)
{
    auto buffer = BridgesBuffers<GRAPH>(graph);
    findBridges(graph, subgraph_mask, starting_vertex, is_bridge, buffer);
}

/// \brief Find, by Tarjan's algorithm, bridges (cut edges) in a subgraph of an undirected graph containing a vertex.
///
/// Tarjan, R. (1974). A note on finding the bridges of a graph.
/// Information Processing Letters 2(6):160-161.
///
/// Runtime complexity O(|E| + |V|).
/// 
/// \param graph An undirected graph.
/// \param subgraph_mask Mask defining the subgraph.
/// \param starting_vertex A vertex inside the subgraph.
/// \param is_bridge Array storing 1 for each edge-index if it is a bridge, or 0 otherwise.
/// \param buffer A pre-allocated buffer object.
///
template<typename GRAPH, typename SUBGRAPH, typename ARRAY>
inline void
findBridges(
    const GRAPH& graph,
    const SUBGRAPH& subgraph_mask,
    std::size_t starting_vertex,
    ARRAY& is_bridge,
    BridgesBuffers<GRAPH>& buffer
)
{
    std::fill(buffer.visited_.begin(), buffer.visited_.end(), 0);

    std::stack<std::size_t> S;

    S.push(starting_vertex);
    buffer.depth_[starting_vertex] = 0;
    buffer.parent_[starting_vertex] = -1;

    while (!S.empty())
    {
        auto v = S.top();
        S.pop();

        if (!buffer.visited_[v])
        {
            buffer.visited_[v] = 1;
            buffer.next_out_arc_[v] = graph.verticesFromVertexBegin(v);
            buffer.min_successor_depth_[v] = buffer.depth_[v];
        }
        else
        {
            auto to = *buffer.next_out_arc_[v];

            if (buffer.min_successor_depth_[to] > buffer.depth_[v])
            {
                typename GRAPH::EdgeIterator e = graph.edgesFromVertexBegin(v) + (buffer.next_out_arc_[v] - graph.verticesFromVertexBegin(v));
                is_bridge[*e] = 1;
            }

            buffer.min_successor_depth_[v] = std::min(buffer.min_successor_depth_[v], buffer.min_successor_depth_[to]);
            ++buffer.next_out_arc_[v];
        }

        while (buffer.next_out_arc_[v] != graph.verticesFromVertexEnd(v))
        {
            typename GRAPH::EdgeIterator e = graph.edgesFromVertexBegin(v) + (buffer.next_out_arc_[v] - graph.verticesFromVertexBegin(v));

            if (
                !subgraph_mask.vertex(*buffer.next_out_arc_[v]) ||
                !subgraph_mask.edge(*e)
                )
            {
                ++buffer.next_out_arc_[v];
                continue;
            }

            auto to = *buffer.next_out_arc_[v];
            if (buffer.visited_[to])
            {
                if(buffer.parent_[v] != to)
                    buffer.min_successor_depth_[v] = std::min(buffer.min_successor_depth_[v], buffer.depth_[to]);

                ++buffer.next_out_arc_[v];
            }
            else
            {
                S.push(v);
                S.push(to);
                buffer.parent_[to] = v;
                buffer.depth_[to] = buffer.depth_[v] + 1;
                is_bridge[*e] = 0;
                break;
            }
        }
    }
}

}
}

#endif
