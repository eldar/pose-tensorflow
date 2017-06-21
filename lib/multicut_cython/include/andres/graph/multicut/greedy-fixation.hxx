#pragma once

#include <cstddef>
#include <iterator>
#include <vector>
#include <algorithm>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <queue>

#include "andres/partition.hxx"

namespace andres {
namespace graph {
namespace multicut {

/// Greedy agglomerative decomposition of a graph.
///
template<typename GRAPH, typename EVA, typename ELA>
void greedyFixation(
    const GRAPH& graph,
    EVA const& edge_values,
    ELA& edge_labels
)
{
    class DynamicGraph
    {
    public:
        DynamicGraph(size_t n) :
            vertices_(n), cut_edges_(n)
        {}

        bool edgeExists(size_t a, size_t b) const
        {
            return !vertices_[a].empty() && vertices_[a].find(b) != vertices_[a].end();
        }

        std::unordered_map<size_t, typename EVA::value_type> const& getAdjacentVertices(size_t v) const
        {
            return vertices_[v];
        }

        typename EVA::value_type getEdgeWeight(size_t a, size_t b) const
        {
            return vertices_[a].at(b);
        }

        bool isCutEdge(size_t a, size_t b) const
        {
            return !cut_edges_[a].empty() && cut_edges_[a].find(b) != cut_edges_[a].end();
        }

        void markEdgeCut(size_t a, size_t b)
        {
            cut_edges_[a].insert(b);
            cut_edges_[b].insert(a);
        }

        void removeVertex(size_t v)
        {
            for (auto& p : vertices_[v])
            {
                vertices_[p.first].erase(v);
                cut_edges_[p.first].erase(v);
            }

            vertices_[v].clear();
            cut_edges_[v].clear();
        }

        void setEdgeWeight(size_t a, size_t b, typename EVA::value_type w)
        {
            vertices_[a][b] = w;
            vertices_[b][a] = w;
        }

    private:
        std::vector<std::unordered_set<size_t>> cut_edges_;
        std::vector<std::unordered_map<size_t, typename EVA::value_type>> vertices_;
    };

    struct Edge
    {
        Edge(size_t _a, size_t _b, typename EVA::value_type _w)
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
        typename EVA::value_type w;

        bool operator <(Edge const& other) const
        {
            return fabs(w) < fabs(other.w);
        }
    };

    std::vector<std::unordered_map<size_t, size_t>> edge_editions(graph.numberOfVertices());
    (graph.numberOfVertices());

    DynamicGraph original_graph_cp(graph.numberOfVertices());
    std::priority_queue<Edge> Q;

    for (size_t i = 0; i < graph.numberOfEdges(); ++i)
    {
        auto a = graph.vertexOfEdge(i, 0);
        auto b = graph.vertexOfEdge(i, 1);

        original_graph_cp.setEdgeWeight(a, b, edge_values[i]);

        auto e = Edge(a, b, edge_values[i]);
        e.edition = ++edge_editions[e.a][e.b];

        Q.push(e);
    }

    andres::Partition<size_t> partition(graph.numberOfVertices());

    while (!Q.empty())
    {
        auto edge = Q.top();
        Q.pop();

        if (!original_graph_cp.edgeExists(edge.a, edge.b) || edge.edition < edge_editions[edge.a][edge.b])
            continue;
        
        if (edge.w > typename EVA::value_type() && !original_graph_cp.isCutEdge(edge.a, edge.b))
        {
            auto stable_vertex = edge.a;
            auto merge_vertex = edge.b;

            if (original_graph_cp.getAdjacentVertices(stable_vertex).size() < original_graph_cp.getAdjacentVertices(merge_vertex).size())
                std::swap(stable_vertex, merge_vertex);

            partition.merge(stable_vertex, merge_vertex);

            for (auto& p : original_graph_cp.getAdjacentVertices(merge_vertex))
            {
                if (p.first == stable_vertex)
                    continue;

                typename EVA::value_type t = typename EVA::value_type();

                if (original_graph_cp.edgeExists(stable_vertex, p.first))
                    t = original_graph_cp.getEdgeWeight(stable_vertex, p.first);

                if (original_graph_cp.isCutEdge(merge_vertex, p.first))
                    original_graph_cp.markEdgeCut(stable_vertex, p.first);

                original_graph_cp.setEdgeWeight(stable_vertex, p.first, p.second + t);

                auto e = Edge(stable_vertex, p.first, p.second + t);
                e.edition = ++edge_editions[e.a][e.b];

                Q.push(e);
            }

            original_graph_cp.removeVertex(merge_vertex);
        }
        else if (edge.w < typename EVA::value_type())
            original_graph_cp.markEdgeCut(edge.a, edge.b);
    }

    for (size_t i = 0; i < graph.numberOfEdges(); ++i)
        edge_labels[i] = partition.find(graph.vertexOfEdge(i, 0)) == partition.find(graph.vertexOfEdge(i, 1)) ? 0 : 1;
}

} // namespace multicut
} // namespace graph
} // namespace andres
