#pragma once
#ifndef ANDRES_GRAPH_MINIMUM_SPANNING_TREE_HXX
#define ANDRES_GRAPH_MINIMUM_SPANNING_TREE_HXX

#include <queue>
#include <stdexcept>
#include <vector>

#include "andres/functional.hxx"
#include "subgraph.hxx"

namespace andres {
namespace graph {

template<typename GRAPH, typename ECA, typename PRED, typename FUNC = Identity<typename ECA::value_type>>
inline typename ECA::value_type 
findMSTPrim(
    const GRAPH&, 
    const ECA&, 
    PRED&, 
    const FUNC& = FUNC()
);

template<typename GRAPH, typename ECA, typename SUBGRAPH, typename PRED, typename FUNC = Identity<typename ECA::value_type>>
inline typename ECA::value_type 
findMSTPrim(
    const GRAPH&, 
    const ECA&, 
    const SUBGRAPH&, 
    PRED&, 
    const FUNC& = FUNC()
);

template<typename GRAPH, typename ECA, typename SUBGRAPH, typename PRED, typename FUNC = Identity<typename ECA::value_type>>
inline typename ECA::value_type 
findMSTPrim(
    const GRAPH&, 
    const ECA&, 
    const SUBGRAPH&, 
    std::size_t, 
    PRED&, 
    const FUNC& = FUNC()
);

template<typename GRAPH, typename ECA, typename PRED, typename FUNC = Identity<typename ECA::value_type>>
inline typename ECA::value_type 
findMSTDynamicProgramming(
    const GRAPH&, 
    const ECA&, 
    PRED&, 
    const FUNC& = FUNC()
);

template<typename GRAPH, typename ECA, typename SUBGRAPH, typename PRED, typename FUNC = Identity<typename ECA::value_type>>
inline typename ECA::value_type 
findMSTDynamicProgramming(
    const GRAPH&, 
    const ECA&, 
    const SUBGRAPH&, 
    PRED&, 
    const FUNC& = FUNC()
);

template<typename GRAPH, typename ECA, typename SUBGRAPH, typename PRED, typename FUNC = Identity<typename ECA::value_type>>
inline typename ECA::value_type 
findMSTDynamicProgramming(
    const GRAPH&, 
    const ECA&, 
    const SUBGRAPH&, 
    std::size_t, 
    PRED&, 
    const FUNC& = FUNC()
);

/// \brief Find, by Prim's algorithm, the minimum spanning forest of an undirected graph.
///
/// Robert C. Prim. (1957). Shortest connection networks and some generalizations.
/// Bell System Technical Journal, 36:1389-1401.
///
/// Runtime complexity O(|E|log|V|).
/// 
/// \param graph An undirected graph.
/// \param edge_weights Edge weights.
/// \param predecessor Vector storing, for each vertex, the index of the edge connecting this vertex to the MST, or graph.numberOfEdges(), for the root of the MST.
/// \param f Functor transforming edge weights.
///
template<typename GRAPH, typename ECA, typename PRED, typename FUNC>
inline typename ECA::value_type 
findMSTPrim(
    const GRAPH& graph, 
    const ECA& edge_weights, 
    PRED& predecessor, 
    const FUNC& f
) 
{
    return findMSTPrim(graph, edge_weights, DefaultSubgraphMask<>(), predecessor, f);
}

/// \brief Find, by Prim's algorithm, the minimum spanning forest of a subgraph of an undirected graph.
///
/// Robert C. Prim. (1957). Shortest connection networks and some generalizations.
/// Bell System Technical Journal, 36:1389-1401.
///
/// Runtime complexity O(|E|log|V|).
///
/// \param graph An undirected graph.
/// \param edge_weights Edge weights.
/// \param subgraph_mask Mask defining the subgraph.
/// \param predecessor Vector storing, for each vertex, the index of the edge connecting this vertex to the MST, or graph.numberOfEdges(), for the root of the MST.
/// \param f Functor transforming edge weights.
///
template<typename GRAPH, typename ECA, typename SUBGRAPH, typename PRED, typename FUNC>
inline typename ECA::value_type 
findMSTPrim(
    const GRAPH& graph, 
    const ECA& edge_weights, 
    const SUBGRAPH& subgraph_mask, 
    PRED& predecessor, 
    const FUNC& f
) 
{
    typedef typename ECA::value_type value_type;

    std::fill(std::begin(predecessor), std::end(predecessor), graph.numberOfEdges());

    value_type mst_value = value_type();
    for (std::size_t i = 0; i < graph.numberOfVertices(); ++i)
        if (predecessor[i] == graph.numberOfEdges() && subgraph_mask.vertex(i))
            mst_value += findMSTPrim(graph, edge_weights, subgraph_mask, i, predecessor, f);

    return mst_value;
}

/// \brief Find, by Prim's algorithm, the MST of the component of a subgraph of an undirected graph containing a (root) vertex.
///
/// Robert C. Prim. (1957). Shortest connection networks and some generalizations.
/// Bell System Technical Journal, 36:1389-1401.
///
/// Runtime complexity O(|E|log|V|).
///
/// \param graph An undirected graph.
/// \param edge_weights Edge weights.
/// \param subgraph_mask Mask defining the subgraph.
/// \param starting_vertex Root vertex of the MST. 
/// \param predecessor Vector storing, for each vertex, the index of the edge connecting this vertex to the MST, or graph.numberOfEdges(), for the root of the MST.
/// \param f Functor transforming edge weights.
///
template<typename GRAPH, typename ECA, typename SUBGRAPH, typename PRED, typename FUNC>
inline typename ECA::value_type 
findMSTPrim(
    const GRAPH& graph, 
    const ECA& edge_weights, 
    const SUBGRAPH& subgraph_mask, 
    std::size_t starting_vertex, 
    PRED& predecessor, 
    const FUNC& f
)
{
    typedef typename ECA::value_type value_type;

    struct element
    {
        element(value_type weight, std::size_t vertex_index) :
            weight_(weight), vertex_index_(vertex_index)
        {}

        bool operator<(const element& other) const
        {
            // as std::priority_queue is a max heap, invert comparison's logic
            return weight_ > other.weight_;
        }

        value_type weight_;
        std::size_t vertex_index_;
    };

    std::vector<value_type> min_edge(graph.numberOfVertices(), std::numeric_limits<value_type>::max());
    std::vector<char> visited(graph.numberOfVertices());

    std::priority_queue<element> Q;
    Q.push(element(value_type(), starting_vertex));

    predecessor[starting_vertex] = graph.numberOfEdges();

    value_type mst_value = value_type();

    while (!Q.empty())
    {
        auto v = Q.top();
        Q.pop();

        if (visited[v.vertex_index_])
            continue;

        visited[v.vertex_index_] = 1;
        mst_value += v.weight_;

        auto e = graph.edgesFromVertexBegin(v.vertex_index_);
        for (auto w = graph.verticesFromVertexBegin(v.vertex_index_); w != graph.verticesFromVertexEnd(v.vertex_index_); ++w, ++e)
            if (subgraph_mask.vertex(*w) &&
                subgraph_mask.edge(*e) &&
                !visited[*w] &&
                *w != v.vertex_index_ &&
                f(edge_weights[*e]) < min_edge[*w]
                )
            {
                min_edge[*w] = f(edge_weights[*e]);
                predecessor[*w] = *e;
                Q.push(element(f(edge_weights[*e]), *w));
            }
    }

    return mst_value;
}

/// \brief Find, by dynamic programming, the minimum spanning forest of an undirected graph.
///
/// Runtime complexity O(|V|^2).
///
/// \param graph An undirected graph.
/// \param edge_weights Edge weights.
/// \param predecessor Vector storing, for each vertex, the index of the edge connecting this vertex to the MST, or graph.numberOfEdges(), for the root of the MST.
/// \param f Functor transforming edge weights.
///
template<typename GRAPH, typename ECA, typename PRED, typename FUNC>
inline typename ECA::value_type 
findMSTDynamicProgramming(
    const GRAPH& graph, 
    const ECA& edge_weights, 
    PRED& predecessor, 
    const FUNC& f
)
{
    return findMSTDynamicProgramming(graph, edge_weights, DefaultSubgraphMask<>(), predecessor, f);
}

/// \brief Find, by Prim's algorithm, the minimum spanning forest of a subgraph of an undirected graph.
///
/// Runtime complexity O(|V|^2).
///
/// \param graph An undirected graph.
/// \param edge_weights Edge weights.
/// \param subgraph_mask Mask defining the subgraph.
/// \param predecessor Vector storing, for each vertex, the index of the edge connecting this vertex to the MST, or graph.numberOfEdges(), for the root of the MST.
/// \param f Functor transforming edge weights.
///
template<typename GRAPH, typename ECA, typename SUBGRAPH, typename PRED, typename FUNC>
inline typename ECA::value_type 
findMSTDynamicProgramming(
    const GRAPH& graph, 
    const ECA& edge_weights, 
    const SUBGRAPH& subgraph_mask, 
    PRED& predecessor, 
    const FUNC& f
)
{
    typedef typename ECA::value_type value_type;

    std::fill(std::begin(predecessor), std::end(predecessor), graph.numberOfEdges());

    value_type mst_value = value_type();
    for (std::size_t i = 0; i < graph.numberOfVertices(); ++i)
        if (predecessor[i] == graph.numberOfEdges() && subgraph_mask.vertex(i))
            mst_value += findMSTDynamicProgramming(graph, edge_weights, subgraph_mask, i, predecessor, f);

    return mst_value;
}

/// \brief Find, by dynamic programming, the MST of the component of a subgraph of an undirected graph containing a (root) vertex.
///
/// Runtime complexity O(|V|^2).
///
/// \param graph An undirected graph.
/// \param edge_weights Edge weights.
/// \param subgraph_mask Mask defining the subgraph.
/// \param starting_vertex Root vertex of the MST. 
/// \param predecessor Vector storing, for each vertex, the index of the edge connecting this vertex to the MST, or graph.numberOfEdges(), for the root of the MST.
/// \param f Functor transforming edge weights.
///
template<typename GRAPH, typename ECA, typename SUBGRAPH, typename PRED, typename FUNC>
inline typename ECA::value_type 
findMSTDynamicProgramming(
    const GRAPH& graph, 
    const ECA& edge_weights, 
    const SUBGRAPH& subgraph_mask, 
    std::size_t starting_vertex, 
    PRED& predecessor, 
    const FUNC& f
)
{
    typedef typename ECA::value_type value_type;

    std::vector<value_type> min_edge(graph.numberOfVertices(), std::numeric_limits<value_type>::max());
    std::vector<char> visited(graph.numberOfVertices());

    min_edge[starting_vertex] = value_type();
    predecessor[starting_vertex] = graph.numberOfEdges();

    value_type mst_value = value_type();

    for (std::size_t i = 0; i < graph.numberOfVertices(); ++i)
    {
        int v = -1;

        for (std::size_t j = 0; j < graph.numberOfVertices(); ++j)
            if (subgraph_mask.vertex(j) &&
                !visited[j] &&
                (v == -1 || min_edge[j] < min_edge[v])
                )
                v = j;

        if (v == -1)
            return mst_value;

        mst_value += min_edge[v];
        visited[v] = 1;

        auto e = graph.edgesFromVertexBegin(v);
        for (auto w = graph.verticesFromVertexBegin(v); w != graph.verticesFromVertexEnd(v); ++w, ++e)
            if (subgraph_mask.vertex(*w) &&
                subgraph_mask.edge(*e) &&
                !visited[*w] &&
                *w != v &&
                f(edge_weights[*e]) < min_edge[*w]
                )
            {
                min_edge[*w] = f(edge_weights[*e]);
                predecessor[*w] = *e;
            }
    }

    return mst_value;
}

} // namespace graph
} // namespace andres

#endif

