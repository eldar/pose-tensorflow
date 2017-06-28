#pragma once
#ifndef ANDRES_GRAPH_BFS_HXX
#define ANDRES_GRAPH_BFS_HXX

#include <cassert>
#include <cstddef>
#include <limits>
#include <queue>
#include <type_traits>
#include <vector>

#include "subgraph.hxx"

namespace andres {
namespace graph {

template<typename S = std::size_t>
class BreadthFirstSearchData {
public:
    typedef S size_type;

    static const size_type NOT_VISITED;

    BreadthFirstSearchData(const size_type size)
        :   depth_(size, NOT_VISITED)
        {}
    template<typename GRAPH>
    BreadthFirstSearchData(const GRAPH& graph)
        :   depth_(graph.numberOfVertices(), NOT_VISITED)
        {}
    size_type add(const size_type v, const size_type depth)
        { depth_[v] = depth; queue_.push(v); }
    void clearQueue()
        { queue_ = std::queue<size_type>(); }
    bool empty() const
        { return queue_.empty(); }
    void markAllNotvisited()
        { std::fill(depth_.begin(), depth_.end(), NOT_VISITED); }
    size_type next()
        { const size_type v = queue_.front(); queue_.pop(); return v; }
    size_type depth(const size_type v) const
        { return depth_[v]; }
    size_type& depth(const size_type v)
        { return depth_[v]; }

private:
    std::vector<size_type> depth_;
    std::queue<size_type> queue_;
};
template<typename S>
   const S BreadthFirstSearchData<S>::NOT_VISITED = std::numeric_limits<S>::max();

template<typename GRAPH, typename CALLBACK>
inline void
breadthFirstSearch(
    const GRAPH& g,
    const std::size_t start_vertex,
    CALLBACK& callback
)
{
    BreadthFirstSearchData<std::size_t> data(g);
    breadthFirstSearch(g, DefaultSubgraphMask<>(), start_vertex, callback, data);
}

template<typename GRAPH, typename CALLBACK, typename = typename std::enable_if<!std::is_lvalue_reference<CALLBACK>::value>::type>
inline void
breadthFirstSearch(
    const GRAPH& g,
    const std::size_t start_vertex,
    CALLBACK&& callback
)
{
    BreadthFirstSearchData<std::size_t> data(g);
    breadthFirstSearch(g, DefaultSubgraphMask<>(), start_vertex, callback, data);
}

template<typename GRAPH, typename CALLBACK>
inline void
breadthFirstSearch(
    const GRAPH& g,
    const std::size_t start_vertex,
    CALLBACK& callback,
    BreadthFirstSearchData<std::size_t>& data
)
{
    breadthFirstSearch(g, DefaultSubgraphMask<>(), start_vertex, callback, data);
}

template<typename GRAPH, typename CALLBACK, typename = typename std::enable_if<!std::is_lvalue_reference<CALLBACK>::value>::type>
inline void
breadthFirstSearch(
    const GRAPH& g,
    const std::size_t start_vertex,
    CALLBACK&& callback,
    BreadthFirstSearchData<std::size_t>& data
)
{
    breadthFirstSearch(g, DefaultSubgraphMask<>(), start_vertex, callback, data);
}

template<typename GRAPH, typename SUBGRAPH, typename CALLBACK>
inline void
breadthFirstSearch(
    const GRAPH& g,
    const SUBGRAPH& subgraph_mask,
    const std::size_t start_vertex,
    CALLBACK& callback
)
{
    BreadthFirstSearchData<std::size_t> data(g);
    breadthFirstSearch(g, subgraph_mask, start_vertex, callback, data);
}

template<typename GRAPH, typename SUBGRAPH, typename CALLBACK, typename = typename std::enable_if<!std::is_lvalue_reference<CALLBACK>::value>::type>
inline void
breadthFirstSearch(
    const GRAPH& g,
    const SUBGRAPH& subgraph_mask,
    const std::size_t start_vertex,
    CALLBACK&& callback
)
{
    BreadthFirstSearchData<std::size_t> data(g);
    breadthFirstSearch(g, subgraph_mask, start_vertex, callback, data);
}

template<typename GRAPH, typename SUBGRAPH, typename CALLBACK>
inline void
breadthFirstSearch(
    const GRAPH& g,
    const SUBGRAPH& subgraph_mask,
    const std::size_t start_vertex,
    CALLBACK& callback,
    BreadthFirstSearchData<std::size_t>& data
)
{
    assert(start_vertex < g.numberOfVertices());

    {
        bool proceed;
        bool add;
        callback(start_vertex, 0, proceed, add);
        if(!proceed) {
            return;
        }
        if(add) {
            data.add(start_vertex, 0);
        }
    }
    while(!data.empty()) {
        const auto v = data.next();
        const auto depth = data.depth(v) + 1;

        auto e_it = g.edgesFromVertexBegin(v);
        for(auto it = g.verticesFromVertexBegin(v); it != g.verticesFromVertexEnd(v); ++it, ++e_it)
            if(data.depth(*it) == BreadthFirstSearchData<std::size_t>::NOT_VISITED &&
                subgraph_mask.vertex(*it) &&
                subgraph_mask.edge(*e_it))
            {
                bool proceed;
                bool add;

                callback(*it, depth, proceed, add);
                
                if(!proceed)
                {
                    data.clearQueue();
                    return;
                }
                
                if(add)
                    data.add(*it, depth);
            }
    }
}

template<typename GRAPH, typename SUBGRAPH, typename CALLBACK, typename = typename std::enable_if<!std::is_lvalue_reference<CALLBACK>::value>::type>
inline void
breadthFirstSearch(
    const GRAPH& g,
    const SUBGRAPH& subgraph_mask,
    const std::size_t start_vertex,
    CALLBACK&& callback,
    BreadthFirstSearchData<std::size_t>& data
)
{
    breadthFirstSearch(g, subgraph_mask, start_vertex, callback, data);
}

} // namespace graph
} // namespace andres

#endif
