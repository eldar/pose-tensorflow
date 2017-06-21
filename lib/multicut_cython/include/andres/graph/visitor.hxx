#pragma once
#ifndef ANDRES_GRAPH_VISITOR_HXX
#define ANDRES_GRAPH_VISITOR_HXX

#include <cstddef>
#include <iostream>

namespace andres {
namespace graph {

/// Visitors can be used to follow the indices of vertices and edges.
///
/// These indices change due to the insetion and removal of vertices and edges.
///
template<class S = std::size_t>
struct IdleGraphVisitor {
    typedef S size_type;

    IdleGraphVisitor() {}
    void insertVertex(const size_type a) const {}
    void insertVertices(const size_type a, const size_type n) const {}
    void eraseVertex(const size_type a) const {}
    void relabelVertex(const size_type a, const size_type b) const {}
    void insertEdge(const size_type a) const {}
    void eraseEdge(const size_type a) const {}
    void relabelEdge(const size_type a, const size_type b) const {}
};

/// Visitors can be used to follow the indices of vertices and edges.
///
/// These indices change due to the insetion and removal of vertices and edges.
///
template<class S = std::size_t>
struct VerboseGraphVisitor {
    typedef S size_type;

    VerboseGraphVisitor() {}
    void insertVertex(const size_type a) const
        { std::cout << "inserting vertex " << a << std::endl; }
    void insertVertices(const size_type a, const size_type n) const
        { std::cout << "inserting " << n << " vertices, starting from index " << a << std::endl; }
    void eraseVertex(const size_type a) const
        { std::cout << "removing vertex " << a << std::endl; }
    void relabelVertex(const size_type a, const size_type b) const
        { std::cout << "relabeling vertex " << a << ". new label is " << b << std::endl; }
    void insertEdge(const size_type a) const
        { std::cout << "inserting edge " << a << std::endl; }
    void eraseEdge(const size_type a) const
        { std::cout << "removing edge " << a << std::endl; }
    void relabelEdge(const size_type a, const size_type b) const
        { std::cout << "relabeling edge " << a << ". new label is " << b << std::endl; }
};

} // namespace graph
} // namespace andres

#endif // #ifndef ANDRES_GRAPH_VISITOR_HXX
