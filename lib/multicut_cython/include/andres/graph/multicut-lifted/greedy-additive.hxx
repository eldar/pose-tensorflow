#pragma once
#ifndef ANDRES_GRAPH_MULTICUT_LIFTED_GREEDY_ADDITIVE_HXX
#define ANDRES_GRAPH_MULTICUT_LIFTED_GREEDY_ADDITIVE_HXX

#include <cstddef>
#include <iterator>
#include <vector>
#include <algorithm>
#include <map>
#include <queue>

#include "andres/partition.hxx"

namespace andres {
namespace graph {
namespace multicut_lifted {

/// Graph decomposition by greedy additive edge contraction.
template<typename ORIGGRAPH, typename LIFTGRAPH, typename ECA>
std::vector<size_t> greedyAdditiveEdgeContraction(
    ORIGGRAPH const& original_graph,
    LIFTGRAPH const& lifted_graph,
    ECA const& edge_costs,
    size_t number_of_clusters_lower_bound = 1,
    size_t number_of_clusters_upper_bound = std::numeric_limits<size_t>::max()
)
{
    class DynamicGraph
    {
    public:
        DynamicGraph(size_t n) :
            vertices_(n)
        {}

        bool edgeExists(size_t a, size_t b) const
        {
            return !vertices_[a].empty() && vertices_[a].find(b) != vertices_[a].end();
        }

        std::map<size_t, typename ECA::value_type> const& getAdjacentVertices(size_t v) const
        {
            return vertices_[v];
        }

        typename ECA::value_type getEdgeWeight(size_t a, size_t b) const
        {
            return vertices_[a].at(b);
        }

        void removeVertex(size_t v)
        {
            for (auto& p : vertices_[v])
                vertices_[p.first].erase(v);

            vertices_[v].clear();
        }

        void setEdgeWeight(size_t a, size_t b, typename ECA::value_type w)
        {
            vertices_[a][b] = w;
            vertices_[b][a] = w;
        }

    private:
        std::vector<std::map<size_t, typename ECA::value_type>> vertices_;
    };

    struct Edge
    {
        Edge(size_t _a, size_t _b, typename ECA::value_type _w)
        {
            if (_a > _b)
                std::swap(_a, _b);

            a = _a;
            b = _b;

            w = _w;
        }

        size_t a;
        size_t b;
        size_t edition;
        typename ECA::value_type w;

        bool operator <(Edge const& other) const
        {
            return w < other.w;
        }
    };

    std::vector<std::map<size_t, size_t>> edge_editions(original_graph.numberOfVertices());
    DynamicGraph original_graph_cp(original_graph.numberOfVertices());
    DynamicGraph lifted_graph_cp(original_graph.numberOfVertices());
    std::priority_queue<Edge> Q;

    for (size_t i = 0; i < original_graph.numberOfEdges(); ++i)
    {
        auto const a = original_graph.vertexOfEdge(i, 0);
        auto const b = original_graph.vertexOfEdge(i, 1);
        
        original_graph_cp.setEdgeWeight(a, b, 1.);
    }

    for (size_t i = 0; i < lifted_graph.numberOfEdges(); ++i)
    {
        auto const a = lifted_graph.vertexOfEdge(i, 0);
        auto const b = lifted_graph.vertexOfEdge(i, 1);

        lifted_graph_cp.setEdgeWeight(a, b, edge_costs[i]);
        
        if (original_graph_cp.edgeExists(a, b))
        {
            auto e = Edge(a, b, edge_costs[i]);
            e.edition = ++edge_editions[e.a][e.b];

            Q.push(e);
        }
    }

    andres::Partition<size_t> partition(original_graph.numberOfVertices());
    while (!Q.empty() && partition.numberOfSets() > number_of_clusters_lower_bound)
    {
        auto edge = Q.top();
        Q.pop();

        if (!original_graph_cp.edgeExists(edge.a, edge.b) || edge.edition < edge_editions[edge.a][edge.b])
            continue;
        
        if (partition.numberOfSets() <= number_of_clusters_upper_bound && edge.w < typename ECA::value_type())
            break;

        auto stable_vertex = edge.a;
        auto merge_vertex = edge.b;

        if (lifted_graph_cp.getAdjacentVertices(stable_vertex).size() < lifted_graph_cp.getAdjacentVertices(merge_vertex).size())
            std::swap(stable_vertex, merge_vertex);

        partition.merge(stable_vertex, merge_vertex);

        for (auto& p : original_graph_cp.getAdjacentVertices(merge_vertex))
        {
            if (p.first == stable_vertex)
                continue;

            original_graph_cp.setEdgeWeight(stable_vertex, p.first, 1.);
        }

        original_graph_cp.removeVertex(merge_vertex);

        for (auto& p : lifted_graph_cp.getAdjacentVertices(merge_vertex))
        {
            if (p.first == stable_vertex)
                continue;

            auto t = typename ECA::value_type();

            if (lifted_graph_cp.edgeExists(stable_vertex, p.first))
                t = lifted_graph_cp.getEdgeWeight(stable_vertex, p.first);

            lifted_graph_cp.setEdgeWeight(stable_vertex, p.first, p.second + t);
            
            if (original_graph_cp.edgeExists(stable_vertex, p.first))
            {
                auto e = Edge(stable_vertex, p.first, p.second + t);
                e.edition = ++edge_editions[e.a][e.b];

                Q.push(e);
            }
        }

        lifted_graph_cp.removeVertex(merge_vertex);
    }

    std::map<size_t, size_t> indices;
    for (size_t v = 0; v < original_graph.numberOfVertices(); ++v)
        indices[partition.find(v)] = 0;

    size_t cnt = 0;
    for (auto& p : indices)
        p.second = cnt++;

    std::vector<size_t> vertex_labels(original_graph.numberOfVertices());
    for (size_t v = 0; v < original_graph.numberOfVertices(); ++v)
        vertex_labels[v] = indices[partition.find(v)];

    return vertex_labels;
}



// function for back compatability with old interface that works with edge labels
template<typename ORIGGRAPH, typename LIFTGRAPH, typename ECA, typename ELA>
void greedyAdditiveEdgeContraction(
    ORIGGRAPH const& original_graph,
    LIFTGRAPH const& lifted_graph,
    ECA const& edge_costs,
    ELA& edge_labels
)
{
    auto vertex_labels = greedyAdditiveEdgeContraction(original_graph, lifted_graph, edge_costs);

    for (size_t e = 0; e < lifted_graph.numberOfEdges(); ++e)
    {
        auto const v0 = lifted_graph.vertexOfEdge(e, 0);
        auto const v1 = lifted_graph.vertexOfEdge(e, 1);

        edge_labels[e] = (vertex_labels[v0] != vertex_labels[v1]) ? 1 : 0;
    }
}

} // namespace multicut_lifted 
} // namespace graph
} // namespace andres

#endif // #ifndef ANDRES_GRAPH_MULTICUT_LIFTED_GREEDY_ADDITIVE_HXX
