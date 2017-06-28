#pragma once
#ifndef ANDRES_GRAPH_COMPLETE_GRAPH_HXX
#define ANDRES_GRAPH_COMPLETE_GRAPH_HXX

#include <cassert>
#include <cstddef>
#include <cmath>
#include <iterator>

#include "adjacency.hxx"
#include "visitor.hxx"

namespace andres {
namespace graph {

/// Complete graph.
template<typename VISITOR = IdleGraphVisitor<std::size_t> >
class CompleteGraph {
public:
    typedef VISITOR Visitor;
    typedef Adjacency<> AdjacencyType;
    
    // \cond SUPPRESS_DOXYGEN
    class AdjacencyIterator
    :   public std::iterator <
            std::random_access_iterator_tag,
            const AdjacencyType
        > {
    public:
        typedef CompleteGraph<Visitor> GraphType;
        typedef std::iterator <
                std::random_access_iterator_tag,
                const AdjacencyType
            > Base;
        typedef typename Base::iterator_category iterator_category;
        typedef typename Base::value_type value_type;
        typedef typename Base::difference_type difference_type;
        typedef typename Base::pointer pointer;
        typedef typename Base::reference reference;

        AdjacencyIterator();
        AdjacencyIterator(const GraphType&);
        AdjacencyIterator(const GraphType&, const std::size_t);
        AdjacencyIterator(const GraphType&, const std::size_t, const std::size_t);

        // increment and decrement
        AdjacencyIterator& operator+=(const difference_type);
        AdjacencyIterator& operator-=(const difference_type);
        AdjacencyIterator& operator++(); // prefix
        AdjacencyIterator& operator--(); // prefix
        AdjacencyIterator operator++(int); // postfix
        AdjacencyIterator operator--(int); // postfix
        AdjacencyIterator operator+(const difference_type) const;
        AdjacencyIterator operator-(const difference_type) const;
        difference_type operator-(const AdjacencyIterator&) const;

        // comparison
        bool operator==(const AdjacencyIterator&) const;
        bool operator!=(const AdjacencyIterator&) const;
        bool operator<(const AdjacencyIterator&) const;
        bool operator<=(const AdjacencyIterator&) const;
        bool operator>(const AdjacencyIterator&) const;
        bool operator>=(const AdjacencyIterator&) const;

        // access
        reference operator*();
        pointer operator->();
        reference operator[](const std::size_t);

    protected:
        const GraphType* graph_;
        std::size_t vertex_;
        std::size_t adjacencyIndex_;
        AdjacencyType adjacency_;
    };

    class VertexIterator
    :   public AdjacencyIterator {
    public:
        typedef CompleteGraph<Visitor> GraphType;
        typedef AdjacencyIterator Base;
        typedef typename Base::iterator_category iterator_category;
        typedef typename Base::difference_type difference_type;
        typedef const std::size_t value_type;
        typedef value_type* pointer;
        typedef value_type& reference;

        VertexIterator();
        VertexIterator(const VertexIterator&);
        VertexIterator(const AdjacencyIterator&);
        VertexIterator(const GraphType&);
        VertexIterator(const GraphType&, const std::size_t);
        VertexIterator(const GraphType&, const std::size_t, const std::size_t);

        // access
        value_type operator*() const;
        value_type operator[](const std::size_t) const;
    private:
        pointer operator->() const;
    };

    class EdgeIterator
    :   public AdjacencyIterator {
    public:
        typedef CompleteGraph<Visitor> GraphType;
        typedef AdjacencyIterator Base;
        typedef typename Base::iterator_category iterator_category;
        typedef typename Base::difference_type difference_type;
        typedef const std::size_t value_type;
        typedef value_type* pointer;
        typedef value_type& reference;

        EdgeIterator();
        EdgeIterator(const EdgeIterator&);
        EdgeIterator(const AdjacencyIterator&);
        EdgeIterator(const GraphType&);
        EdgeIterator(const GraphType&, const std::size_t);
        EdgeIterator(const GraphType&, const std::size_t, const std::size_t);

        // access
        value_type operator*() const;
        value_type operator[](const std::size_t) const;
    private:
        pointer operator->() const;
    };
    // \endcond

    // construction
    CompleteGraph(const Visitor& = Visitor());
    CompleteGraph(const std::size_t, const Visitor& = Visitor());
    void assign(const Visitor& = Visitor());
    void assign(const std::size_t, const Visitor& = Visitor());

    // iterator access (compatible with Digraph)
    VertexIterator verticesFromVertexBegin(const std::size_t) const;
    VertexIterator verticesFromVertexEnd(const std::size_t) const;
    VertexIterator verticesToVertexBegin(const std::size_t) const;
    VertexIterator verticesToVertexEnd(const std::size_t) const;
    EdgeIterator edgesFromVertexBegin(const std::size_t) const;
    EdgeIterator edgesFromVertexEnd(const std::size_t) const;
    EdgeIterator edgesToVertexBegin(const std::size_t) const;
    EdgeIterator edgesToVertexEnd(const std::size_t) const;
    AdjacencyIterator adjacenciesFromVertexBegin(const std::size_t) const;
    AdjacencyIterator adjacenciesFromVertexEnd(const std::size_t) const;
    AdjacencyIterator adjacenciesToVertexBegin(const std::size_t) const;
    AdjacencyIterator adjacenciesToVertexEnd(const std::size_t) const;

    // access (compatible with Digraph)
    std::size_t numberOfVertices() const;
    std::size_t numberOfEdges() const;
    std::size_t numberOfEdgesFromVertex(const std::size_t) const;
    std::size_t numberOfEdgesToVertex(const std::size_t) const;
    std::size_t vertexOfEdge(const std::size_t, const std::size_t) const;
    std::size_t edgeFromVertex(const std::size_t, const std::size_t) const;
    std::size_t edgeToVertex(const std::size_t, const std::size_t) const;
    std::size_t vertexFromVertex(const std::size_t, const std::size_t) const;
    std::size_t vertexToVertex(const std::size_t, const std::size_t) const;
    AdjacencyType adjacencyFromVertex(const std::size_t, const std::size_t) const;
    AdjacencyType adjacencyToVertex(const std::size_t, const std::size_t) const;
    std::pair<bool, std::size_t> findEdge(const std::size_t, const std::size_t) const;
    bool multipleEdgesEnabled() const;

private:
    std::size_t edgeOfStrictlyIncreasingPairOfVertices(const std::size_t, const std::size_t) const;

    std::size_t numberOfVertices_;
    Visitor visitor_;
};

/// Construct a complete graph.
///
/// \param visitor Visitor to follow changes of integer indices of vertices and edges.
///
template<typename VISITOR>
inline
CompleteGraph<VISITOR>::CompleteGraph(
    const Visitor& visitor
)
:   numberOfVertices_(0),
    visitor_(visitor)
{}

/// Construct a complete graph with an initial number of vertices.
///
/// \param numberOfVertices Number of vertices.
/// \param visitor Visitor to follow changes of integer indices of vertices and edges.
///
template<typename VISITOR>
inline
CompleteGraph<VISITOR>::CompleteGraph(
    const std::size_t numberOfVertices,
    const Visitor& visitor
)
:   numberOfVertices_(numberOfVertices),
    visitor_(visitor)
{}

/// Clear a complete graph.
///
/// \param visitor Visitor to follow changes of integer indices of vertices and edges.
///
template<typename VISITOR>
inline void
CompleteGraph<VISITOR>::assign(
    const Visitor& visitor
) {
    numberOfVertices_ = 0;
    visitor_ = visitor;
}

/// Clear a complete graph with an initial number of vertices.
///
/// \param numberOfVertices Number of vertices.
/// \param visitor Visitor to follow changes of integer indices of vertices and edges.
///
template<typename VISITOR>
inline void
CompleteGraph<VISITOR>::assign(
    const std::size_t numberOfVertices,
    const Visitor& visitor
) {
    numberOfVertices_ = numberOfVertices;
    visitor_ = visitor;
}

/// Get an iterator to the beginning of the sequence of vertices reachable from a given vertex via a single edge.
///
/// \param vertex Integer index of the vertex.
/// \return VertexIterator.
///
/// \sa verticesFromVertexEnd()
///
template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::VertexIterator
CompleteGraph<VISITOR>::verticesFromVertexBegin(
    const std::size_t vertex
) const {
    return VertexIterator(*this, vertex, 0);
}

/// Get an iterator to the end of the sequence of vertices reachable from a given vertex via a single edge.
///
/// \param vertex Integer index of the vertex.
/// \return VertexIterator.
///
/// \sa verticesFromVertexBegin()
///
template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::VertexIterator
CompleteGraph<VISITOR>::verticesFromVertexEnd(
    const std::size_t vertex
) const {
    return VertexIterator(*this, vertex, numberOfEdgesFromVertex(vertex));
}

/// Get an iterator to the beginning of the sequence of vertices from which a given vertex is reachable via a single edge.
///
/// \param vertex Integer index of the vertex.
/// \return VertexIterator.
///
/// \sa verticesToVertexEnd()
///
template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::VertexIterator
CompleteGraph<VISITOR>::verticesToVertexBegin(
    const std::size_t vertex
) const {
    return VertexIterator(*this, vertex, 0);
}

/// Get an iterator to the end of the sequence of vertices from which a given vertex is reachable via a single edge.
///
/// \param vertex Integer index of the vertex.
/// \return VertexIterator.
///
/// \sa verticesToVertexBegin()
///
template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::VertexIterator
CompleteGraph<VISITOR>::verticesToVertexEnd(
    const std::size_t vertex
) const {
    return VertexIterator(*this, vertex, numberOfEdgesFromVertex(vertex));
}

/// Get an iterator to the beginning of the sequence of edges that originate from a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return EdgeIterator.
///
/// \sa edgesFromVertexEnd()
///
template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::EdgeIterator
CompleteGraph<VISITOR>::edgesFromVertexBegin(
    const std::size_t vertex
) const {
    return EdgeIterator(*this, vertex, 0);
}

/// Get an iterator to the end of the sequence of edges that originate from a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return EdgeIterator.
///
/// \sa edgesFromVertexBegin()
///
template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::EdgeIterator
CompleteGraph<VISITOR>::edgesFromVertexEnd(
    const std::size_t vertex
) const {
    return EdgeIterator(*this, vertex, numberOfEdgesFromVertex(vertex));
}

/// Get an iterator to the beginning of the sequence of edges that are incident to a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return EdgeIterator.
///
/// \sa edgesToVertexEnd()
///
template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::EdgeIterator
CompleteGraph<VISITOR>::edgesToVertexBegin(
    const std::size_t vertex
) const {
    return EdgeIterator(*this, vertex, 0);
}

/// Get an iterator to the end of the sequence of edges that are incident to a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return EdgeIterator.
///
/// \sa edgesToVertexBegin()
///
template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::EdgeIterator
CompleteGraph<VISITOR>::edgesToVertexEnd(
    const std::size_t vertex
) const {
    return EdgeIterator(*this, vertex, numberOfEdgesFromVertex(vertex));
}

/// Get an iterator to the beginning of the sequence of adjacencies that originate from a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return AdjacencyIterator.
///
/// \sa adjacenciesFromVertexEnd()
///
template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::AdjacencyIterator
CompleteGraph<VISITOR>::adjacenciesFromVertexBegin(
    const std::size_t vertex
) const {
    return AdjacencyIterator(*this, vertex, 0);
}

/// Get an iterator to the end of the sequence of adjacencies that originate from a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return AdjacencyIterator.
///
/// \sa adjacenciesFromVertexBegin()
///
template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::AdjacencyIterator
CompleteGraph<VISITOR>::adjacenciesFromVertexEnd(
    const std::size_t vertex
) const {
    return AdjacencyIterator(*this, vertex, numberOfEdgesFromVertex(vertex));
}

/// Get an iterator to the beginning of the sequence of adjacencies incident to a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return AdjacencyIterator.
///
/// \sa adjacenciesToVertexEnd()
///
template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::AdjacencyIterator
CompleteGraph<VISITOR>::adjacenciesToVertexBegin(
    const std::size_t vertex
) const {
    return AdjacencyIterator(*this, vertex, 0);
}

/// Get an iterator to the end of the sequence of adjacencies incident to a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return AdjacencyIterator.
///
/// \sa adjacenciesToVertexBegin()
///
template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::AdjacencyIterator
CompleteGraph<VISITOR>::adjacenciesToVertexEnd(
    const std::size_t vertex
) const {
    return AdjacencyIterator(*this, vertex, numberOfEdgesFromVertex(vertex));
}

/// Get the number of vertices.
///
template<typename VISITOR>
inline std::size_t
CompleteGraph<VISITOR>::numberOfVertices() const {
    return numberOfVertices_;
}

/// Get the number of edges.
///
template<typename VISITOR>
inline std::size_t
CompleteGraph<VISITOR>::numberOfEdges() const {
    return numberOfVertices() * (numberOfVertices() - 1) / 2;
}

/// Get the number of edges that originate from a given vertex.
///
/// \param vertex Integer index of a vertex.
///
/// \sa edgeFromVertex()
///
template<typename VISITOR>
inline std::size_t
CompleteGraph<VISITOR>::numberOfEdgesFromVertex(
    const std::size_t vertex
) const {
    assert(vertex < numberOfVertices());
    return numberOfVertices() - 1;
}

/// Get the number of edges that are incident to a given vertex.
///
/// \param vertex Integer index of a vertex.
///
/// \sa edgeToVertex()
///
template<typename VISITOR>
inline std::size_t
CompleteGraph<VISITOR>::numberOfEdgesToVertex(
    const std::size_t vertex
) const {
    assert(vertex < numberOfVertices());
    return numberOfVertices() - 1;
}

/// Get the integer index of a vertex of an edge.
///
/// \param edge Integer index of an edge.
/// \param j Number of the vertex in the edge; either 0 or 1.
///
template<typename VISITOR>
inline std::size_t
CompleteGraph<VISITOR>::vertexOfEdge(
    const std::size_t edge,
    const std::size_t j
) const {
    assert(edge < numberOfEdges());
    assert(j < 2);
    double p = static_cast<double>(numberOfVertices() * 2 - 1) / 2;
    double q = static_cast<double>(edge) * 2;
    std::size_t vertex0 = static_cast<std::size_t>(p - std::sqrt(p * p - q));
    if(j == 0) {
        return vertex0;
    }
    else {
        return edge + vertex0 * (vertex0 + 1) / 2 - (numberOfVertices() - 1) * vertex0 + 1;
    }
}

/// Get the integer index of an edge that originates from a given vertex.
///
/// \param vertex Integer index of a vertex.
/// \param j Number of the edge; between 0 and numberOfEdgesFromVertex(vertex) - 1.
///
/// \sa numberOfEdgesFromVertex()
///
template<typename VISITOR>
inline std::size_t
CompleteGraph<VISITOR>::edgeFromVertex(
    const std::size_t vertex,
    const std::size_t j
) const {
    assert(j < numberOfEdgesFromVertex(vertex));
    if(j < vertex) {
        const std::size_t vertexAdjacent = j;
        const std::size_t edgeAdjacent = edgeOfStrictlyIncreasingPairOfVertices(vertexAdjacent, vertex);
        return edgeAdjacent;
    }
    else {
        const std::size_t vertexAdjacent = j + 1;
        const std::size_t edgeAdjacent = edgeOfStrictlyIncreasingPairOfVertices(vertex, vertexAdjacent);
        return edgeAdjacent;
    }
}

/// Get the integer index of an edge that is incident to a given vertex.
///
/// \param vertex Integer index of a vertex.
/// \param j Number of the edge; between 0 and numberOfEdgesFromVertex(vertex) - 1.
///
/// \sa numberOfEdgesToVertex()
///
template<typename VISITOR>
inline std::size_t
CompleteGraph<VISITOR>::edgeToVertex(
    const std::size_t vertex,
    const std::size_t j
) const {
    assert(j < numberOfEdgesToVertex(vertex));
    return edgeFromVertex(vertex, j);
}

/// Get the integer index of a vertex reachable from a given vertex via a single edge.
///
/// \param vertex Integer index of a vertex.
/// \param j Number of the vertex; between 0 and numberOfEdgesFromVertex(vertex) - 1.
///
/// \sa numberOfEdgesFromVertex()
///
template<typename VISITOR>
inline std::size_t
CompleteGraph<VISITOR>::vertexFromVertex(
    const std::size_t vertex,
    const std::size_t j
) const {
    assert(j < numberOfEdgesFromVertex(vertex));
    if(j < vertex) {
        return j;
    }
    else {
        return j + 1;
    }
}

/// Get the integer index of a vertex from which a given vertex is reachable via a single edge.
///
/// \param vertex Integer index of a vertex.
/// \param j Number of the vertex; between 0 and numberOfEdgesFromVertex(vertex) - 1.
///
/// \sa numberOfEdgesFromVertex()
///
template<typename VISITOR>
inline std::size_t
CompleteGraph<VISITOR>::vertexToVertex(
    const std::size_t vertex,
    const std::size_t j
) const {
    assert(j < numberOfEdgesToVertex(vertex));
    return vertexFromVertex(vertex, j);
}

/// Get the j-th adjacency from a vertex.
///
/// \param vertex Vertex.
/// \param j Number of the adjacency.
///
template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::AdjacencyType
CompleteGraph<VISITOR>::adjacencyFromVertex(
    const std::size_t vertex,
    const std::size_t j
) const {
    assert(j < numberOfEdgesToVertex(vertex));
    if(j < vertex) {
        const std::size_t vertexAdjacent = j;
        const std::size_t edgeAdjacent = edgeOfStrictlyIncreasingPairOfVertices(vertexAdjacent, vertex);
        return AdjacencyType(vertexAdjacent, edgeAdjacent);
    }
    else {
        const std::size_t vertexAdjacent = j + 1;
        const std::size_t edgeAdjacent = edgeOfStrictlyIncreasingPairOfVertices(vertex, vertexAdjacent);
        return AdjacencyType(vertexAdjacent, edgeAdjacent);
    }
}

/// Get the j-th adjacency to a vertex.
///
/// \param vertex Vertex.
/// \param j Number of the adjacency.
///
template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::AdjacencyType
CompleteGraph<VISITOR>::adjacencyToVertex(
    const std::size_t vertex,
    const std::size_t j
) const {
    return adjacencyFromVertex(vertex, j);
}

/// Search for an edge (in constant time).
///
/// Indexing: findEdge(vertex0, vertex1)
///    - 0 1 2
///    0 - 3 4
///    1 3 - 5
///    2 4 5 -
///
/// \param vertex0 first vertex of the edge.
/// \param vertex1 second vertex of the edge.
/// \return if an edge from vertex0 to vertex1 exists, pair.first is true
///     and pair.second is the index of such an edge. if no edge from vertex0
///     to vertex1 exists, pair.first is false and pair.second is undefined.
///
template<typename VISITOR>
inline std::pair<bool, std::size_t>
CompleteGraph<VISITOR>::findEdge(
    const std::size_t vertex0,
    const std::size_t vertex1
) const {
    assert(vertex0 < numberOfVertices());
    assert(vertex1 < numberOfVertices());
    if(vertex0 == vertex1) {
        return std::pair<bool, std::size_t>(false, 0);
    }
    else if(vertex0 < vertex1) {
        return std::pair<bool, std::size_t>(true, edgeOfStrictlyIncreasingPairOfVertices(vertex0, vertex1));
    }
    else {
        return std::pair<bool, std::size_t>(true, edgeOfStrictlyIncreasingPairOfVertices(vertex1, vertex0));
    }
}

/// Indicate if multiple edges are enabled.
///
/// \return false
///
template<typename VISITOR>
inline bool
CompleteGraph<VISITOR>::multipleEdgesEnabled() const {
    return false;
}

// private
template<typename VISITOR>
inline std::size_t
CompleteGraph<VISITOR>::edgeOfStrictlyIncreasingPairOfVertices(
    const std::size_t vertex0,
    const std::size_t vertex1
) const {
    assert(vertex1 < numberOfVertices());
    assert(vertex0 < vertex1);
    return (numberOfVertices() - 1) * vertex0 - vertex0 * (vertex0 + 1) / 2 + vertex1 - 1;
}

// \cond SUPPRESS_DOXYGEN

// implementation of AdjacencyIterator

template<typename VISITOR>
inline
CompleteGraph<VISITOR>::AdjacencyIterator::AdjacencyIterator()
:   graph_(0),
    vertex_(0),
    adjacencyIndex_(0),
    adjacency_()
{}

template<typename VISITOR>
inline
CompleteGraph<VISITOR>::AdjacencyIterator::AdjacencyIterator(
    const GraphType& graph
)
:   graph_(&graph),
    vertex_(0),
    adjacencyIndex_(0),
    adjacency_()
{}

template<typename VISITOR>
inline
CompleteGraph<VISITOR>::AdjacencyIterator::AdjacencyIterator(
    const GraphType& graph,
    const std::size_t vertex
)
:   graph_(&graph),
    vertex_(vertex),
    adjacencyIndex_(0),
    adjacency_()
{
    assert(vertex < graph.numberOfVertices());
}

template<typename VISITOR>
inline
CompleteGraph<VISITOR>::AdjacencyIterator::AdjacencyIterator(
    const GraphType& graph,
    const std::size_t vertex,
    const std::size_t adjacencyIndex
)
:   graph_(&graph),
    vertex_(vertex),
    adjacencyIndex_(adjacencyIndex),
    adjacency_()
{
    assert(vertex < graph.numberOfVertices());
    assert(adjacencyIndex <= graph.numberOfVertices());
}

template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::AdjacencyIterator&
CompleteGraph<VISITOR>::AdjacencyIterator::operator+=(
    const difference_type d
) {
    adjacencyIndex_ += d;
    return *this;
}

template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::AdjacencyIterator&
CompleteGraph<VISITOR>::AdjacencyIterator::operator-=(
    const difference_type d
) {
    adjacencyIndex_ -= d;
    return *this;
}

template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::AdjacencyIterator&
CompleteGraph<VISITOR>::AdjacencyIterator::operator++() {
    ++adjacencyIndex_;
    return *this;
}

template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::AdjacencyIterator&
CompleteGraph<VISITOR>::AdjacencyIterator::operator--() {
    --adjacencyIndex_;
    return *this;
}

template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::AdjacencyIterator
CompleteGraph<VISITOR>::AdjacencyIterator::operator++(int) {
    AdjacencyIterator copy = *this;
    ++adjacencyIndex_;
    return copy;
}

template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::AdjacencyIterator
CompleteGraph<VISITOR>::AdjacencyIterator::operator--(int) {
    AdjacencyIterator copy = *this;
    --adjacencyIndex_;
    return copy;
}

template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::AdjacencyIterator
CompleteGraph<VISITOR>::AdjacencyIterator::operator+(
    const difference_type d
) const {
    AdjacencyIterator copy = *this;
    copy += d;
    return copy;
}

template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::AdjacencyIterator
CompleteGraph<VISITOR>::AdjacencyIterator::operator-(
    const difference_type d
) const {
    AdjacencyIterator copy = *this;
    copy -= d;
    return copy;
}

template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::AdjacencyIterator::difference_type
CompleteGraph<VISITOR>::AdjacencyIterator::operator-(
    const AdjacencyIterator& adjacencyIterator
) const {
    return adjacencyIndex_ - adjacencyIterator.adjacencyIndex_;
}

template<typename VISITOR>
inline bool
CompleteGraph<VISITOR>::AdjacencyIterator::operator==(
    const AdjacencyIterator& other
) const {
    return adjacencyIndex_ == other.adjacencyIndex_
        && vertex_ == other.vertex_
        && graph_ == other.graph_;
}

template<typename VISITOR>
inline bool
CompleteGraph<VISITOR>::AdjacencyIterator::operator!=(
    const AdjacencyIterator& other
) const {
    return adjacencyIndex_ != other.adjacencyIndex_
        || vertex_ != other.vertex_
        || graph_ != other.graph_;
}

template<typename VISITOR>
inline bool
CompleteGraph<VISITOR>::AdjacencyIterator::operator<(
    const AdjacencyIterator& other
) const {
    return adjacencyIndex_ < other.adjacencyIndex_
        && vertex_ == other.vertex_
        && graph_ == other.graph_;
}

template<typename VISITOR>
inline bool
CompleteGraph<VISITOR>::AdjacencyIterator::operator<=(
    const AdjacencyIterator& other
) const {
    return adjacencyIndex_ <= other.adjacencyIndex_
        && vertex_ == other.vertex_
        && graph_ == other.graph_;
}


template<typename VISITOR>
inline bool
CompleteGraph<VISITOR>::AdjacencyIterator::operator>(
    const AdjacencyIterator& other
) const {
    return adjacencyIndex_ > other.adjacencyIndex_
        && vertex_ == other.vertex_
        && graph_ == other.graph_;
}

template<typename VISITOR>
inline bool
CompleteGraph<VISITOR>::AdjacencyIterator::operator>=(
    const AdjacencyIterator& other
) const {
    return adjacencyIndex_ >= other.adjacencyIndex_
        && vertex_ == other.vertex_
        && graph_ == other.graph_;
}

template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::AdjacencyIterator::reference
CompleteGraph<VISITOR>::AdjacencyIterator::operator*() {
    adjacency_ = graph_->adjacencyFromVertex(vertex_, adjacencyIndex_);
    return adjacency_;
}

template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::AdjacencyIterator::pointer
CompleteGraph<VISITOR>::AdjacencyIterator::operator->() {
    adjacency_ = graph_->adjacencyFromVertex(vertex_, adjacencyIndex_);
    return &adjacency_;
}

template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::AdjacencyIterator::reference
CompleteGraph<VISITOR>::AdjacencyIterator::operator[](
    const std::size_t j
) {
    adjacency_ = graph_->adjacencyFromVertex(vertex_, adjacencyIndex_ + j);
    return adjacency_;
}

// implementation of VertexIterator

template<typename VISITOR>
inline
CompleteGraph<VISITOR>::VertexIterator::VertexIterator()
:   Base()
{}

template<typename VISITOR>
inline
CompleteGraph<VISITOR>::VertexIterator::VertexIterator(
    const GraphType& graph
)
:   Base(graph)
{}

template<typename VISITOR>
inline
CompleteGraph<VISITOR>::VertexIterator::VertexIterator(
    const GraphType& graph,
    const std::size_t vertex
)
:   Base(graph, vertex)
{}

template<typename VISITOR>
inline
CompleteGraph<VISITOR>::VertexIterator::VertexIterator(
    const GraphType& graph,
    const std::size_t vertex,
    const std::size_t adjacencyIndex
)
:   Base(graph, vertex, adjacencyIndex)
{}

template<typename VISITOR>
inline
CompleteGraph<VISITOR>::VertexIterator::VertexIterator(
    const VertexIterator& it
)
:   Base(*(it.graph_), it.vertex_, it.adjacencyIndex_)
{}

template<typename VISITOR>
inline
CompleteGraph<VISITOR>::VertexIterator::VertexIterator(
    const AdjacencyIterator& it
)
:   Base(it)
{}

template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::VertexIterator::value_type
CompleteGraph<VISITOR>::VertexIterator::operator*() const {
    return Base::graph_->vertexFromVertex(Base::vertex_, Base::adjacencyIndex_);
}

template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::VertexIterator::value_type
CompleteGraph<VISITOR>::VertexIterator::operator[](
    const std::size_t j
) const {
    return Base::graph_->vertexFromVertex(Base::vertex_, Base::adjacencyIndex_ + j);
}

// implementation of EdgeIterator

template<typename VISITOR>
inline
CompleteGraph<VISITOR>::EdgeIterator::EdgeIterator()
:   Base()
{}

template<typename VISITOR>
inline
CompleteGraph<VISITOR>::EdgeIterator::EdgeIterator(
    const GraphType& graph
)
:   Base(graph)
{}

template<typename VISITOR>
inline
CompleteGraph<VISITOR>::EdgeIterator::EdgeIterator(
    const GraphType& graph,
    const std::size_t vertex
)
:   Base(graph, vertex)
{}

template<typename VISITOR>
inline
CompleteGraph<VISITOR>::EdgeIterator::EdgeIterator(
    const GraphType& graph,
    const std::size_t vertex,
    const std::size_t adjacencyIndex
)
:   Base(graph, vertex, adjacencyIndex)
{}

template<typename VISITOR>
inline
CompleteGraph<VISITOR>::EdgeIterator::EdgeIterator(
    const EdgeIterator& it
)
:   Base(*(it.graph_), it.vertex_, it.adjacencyIndex_)
{}

template<typename VISITOR>
inline
CompleteGraph<VISITOR>::EdgeIterator::EdgeIterator(
    const AdjacencyIterator& it
)
:   Base(it)
{}

template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::EdgeIterator::value_type
CompleteGraph<VISITOR>::EdgeIterator::operator*() const {
    return Base::graph_->edgeFromVertex(Base::vertex_, Base::adjacencyIndex_);
}

template<typename VISITOR>
inline typename CompleteGraph<VISITOR>::EdgeIterator::value_type
CompleteGraph<VISITOR>::EdgeIterator::operator[](
    const std::size_t j
) const {
    return Base::graph_->edgeFromVertex(Base::vertex_, Base::adjacencyIndex_ + j);
}

// \endcond

} // namespace graph
} // namespace andres

#endif // #ifndef ANDRES_GRAPH_COMPLETE_GRAPH_HXX
