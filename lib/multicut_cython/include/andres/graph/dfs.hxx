#pragma once
#ifndef ANDRES_GRAPH_DFS_HXX
#define ANDRES_GRAPH_DFS_HXX

#include <cassert>
#include <cstddef>
#include <stack>
#include <type_traits>
#include <vector>

#include "subgraph.hxx"

namespace andres {
namespace graph {

template<typename S = std::size_t>
class DepthFirstSearchData {
public:
    typedef S size_type;

    DepthFirstSearchData(const size_type size)
        :   visited_(size)
        {}
    template<typename GRAPH>
    DepthFirstSearchData(const GRAPH& graph)
        :   visited_(graph.numberOfVertices())
        {}
    size_type add(const size_type v)
        { stack_.push(v); }
    void clearStack()
        { stack_ = std::stack<size_type>(); }
    bool empty() const
        { return stack_.empty(); }
    void markAllNotvisited()
        { std::fill(visited_.begin(), visited_.end(), 0); }
    size_type next()
        { const size_type v = stack_.top(); stack_.pop(); return v; }
    unsigned char& visited(const size_type v)
        { return visited_[v]; }
    unsigned char visited(const size_type v) const
        { return visited_[v]; }

private:
    std::vector<unsigned char> visited_;
    std::stack<size_type> stack_;
};

template<typename GRAPH, typename CALLBACK>
inline void
depthFirstSearch(
    const GRAPH& g,
    const std::size_t start_vertex,
    CALLBACK& callback
)
{
    DepthFirstSearchData<std::size_t> data(g);
    depthFirstSearch(g, DefaultSubgraphMask<>(), start_vertex, callback, data);
}

template<typename GRAPH, typename CALLBACK, typename = typename std::enable_if<!std::is_lvalue_reference<CALLBACK>::value>::type>
inline void
depthFirstSearch(
    const GRAPH& g,
    const std::size_t start_vertex,
    CALLBACK&& callback
)
{
    DepthFirstSearchData<std::size_t> data(g);
    depthFirstSearch(g, DefaultSubgraphMask<>(), start_vertex, callback, data);
}

template<typename GRAPH, typename CALLBACK>
inline void
depthFirstSearch(
    const GRAPH& g,
    const std::size_t start_vertex,
    CALLBACK& callback,
    DepthFirstSearchData<std::size_t>& data
)
{
    depthFirstSearch(g, DefaultSubgraphMask<>(), start_vertex, callback, data);
}

template<typename GRAPH, typename CALLBACK, typename = typename std::enable_if<!std::is_lvalue_reference<CALLBACK>::value>::type>
inline void
depthFirstSearch(
    const GRAPH& g,
    const std::size_t start_vertex,
    CALLBACK&& callback,
    DepthFirstSearchData<std::size_t>& data
)
{
    depthFirstSearch(g, DefaultSubgraphMask<>(), start_vertex, callback, data);
}

template<typename GRAPH, typename SUBGRAPH, typename CALLBACK>
inline void
depthFirstSearch(
    const GRAPH& g,
    const SUBGRAPH& subgraph_mask,
    const std::size_t start_vertex,
    CALLBACK& callback
)
{
    DepthFirstSearchData<std::size_t> data(g);
    depthFirstSearch(g, subgraph_mask, start_vertex, callback, data);
}

template<typename GRAPH, typename SUBGRAPH, typename CALLBACK, typename = typename std::enable_if<!std::is_lvalue_reference<CALLBACK>::value>::type>
inline void
depthFirstSearch(
    const GRAPH& g,
    const SUBGRAPH& subgraph_mask,
    const std::size_t start_vertex,
    CALLBACK&& callback
)
{
    DepthFirstSearchData<std::size_t> data(g);
    depthFirstSearch(g, subgraph_mask, start_vertex, callback, data);
}

template<typename GRAPH, typename SUBGRAPH, typename CALLBACK>
inline void
depthFirstSearch(
    const GRAPH& g,
    const SUBGRAPH& subgraph_mask,
    const std::size_t start_vertex,
    CALLBACK& callback,
    DepthFirstSearchData<std::size_t>& data
)
{
    assert(start_vertex < g.numberOfVertices());

    data.add(start_vertex);

    while(!data.empty())
    {
        auto v = data.next();

        if (!data.visited(v))
        {
            data.visited(v) = 1;
            
            bool proceed;
            bool addNeighbors;

            callback(v, proceed, addNeighbors);
            
            if (!proceed)
            {
                data.clearStack();
                return;
            }
            
            if (addNeighbors)
            {
                auto e_it = g.edgesFromVertexBegin(v);
                for(auto it = g.verticesFromVertexBegin(v); it != g.verticesFromVertexEnd(v); ++it, ++e_it)
                    if (!data.visited(*it) && subgraph_mask.vertex(*it) && subgraph_mask.edge(*e_it))
                        data.add(*it);
            }
        }
    }
}

template<typename GRAPH, typename SUBGRAPH, typename CALLBACK, typename = typename std::enable_if<!std::is_lvalue_reference<CALLBACK>::value>::type>
inline void
depthFirstSearch(
    const GRAPH& g,
    const SUBGRAPH& subgraph_mask,
    const std::size_t start_vertex,
    CALLBACK&& callback,
    DepthFirstSearchData<std::size_t>& data
)
{
    depthFirstSearch(g, subgraph_mask, start_vertex, callback, data);
}

} // namespace graph
} // namespace andres

#endif
