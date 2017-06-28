
#pragma once
#ifndef ANDRES_GRAPH_GRID_GRAPH_HXX
#define ANDRES_GRAPH_GRID_GRAPH_HXX

#include <cassert>
#include <cstddef>
#include <cmath>
#include <iterator>
#include <array>
#include <initializer_list>

#include "adjacency.hxx"
#include "visitor.hxx"

namespace andres {

namespace graph {

/// D-dimensional grid graph.
template <unsigned char D = 2, class VISITOR = IdleGraphVisitor<std::size_t> >
class GridGraph {
public:
    typedef std::size_t size_type;
    typedef VISITOR Visitor;
    typedef andres::graph::Adjacency<size_type> AdjacencyType;

    static const size_type DIMENSION = static_cast<size_type> (D);

    typedef std::array<size_type, DIMENSION> VertexCoordinate;

    /// \class EdgeCoordinate
    /// Describes an edge as the integer index of the minimum of the
    /// two endpoints and the direction along which it is drawn.
    struct EdgeCoordinate {
        EdgeCoordinate(const VertexCoordinate&, const size_type, bool = false);
        EdgeCoordinate();
        /// The minimum of the two endpoints of the edge is specified
        /// as an integer and accessed by the \b pivot member.
        /// This specifies the \e origin of the edge, in a direction
        /// towards the positive orthant.
        /// (i.e.: in specifying the pivot of an edge, an isSmaller
        /// of false is implied.)
        VertexCoordinate pivotCoordinate_;
        /// The dimension along which the edge is drawn.
        size_type dimension_;
    };

    /// AdjacencyIterator
    // \cond SUPPRESS_DOXYGEN
    class AdjacencyIterator
        :   public std::iterator <
                std::random_access_iterator_tag,
                const AdjacencyType
            >  {
    public:        
        typedef GridGraph<DIMENSION, Visitor> GraphType;
        typedef std::iterator <
                std::random_access_iterator_tag,
                const AdjacencyType
            > Base;
        typedef typename Base::difference_type difference_type;
        typedef typename Base::pointer pointer;
        typedef typename Base::reference reference;

        AdjacencyIterator();
        AdjacencyIterator(const GraphType&);
        AdjacencyIterator(const GraphType&, const size_type);
        AdjacencyIterator(const GraphType&, const size_type, const size_type);

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
        reference operator[](const difference_type);

    protected:
        const GraphType* graph_;
        size_type vertex_;
        size_type adjacencyIndex_;
        AdjacencyType adjacency_;
    };

    class VertexIterator
            :   public AdjacencyIterator {
    public:
        typedef GridGraph<DIMENSION, Visitor> GraphType;
        typedef AdjacencyIterator Base;
        typedef const size_type value_type;
        typedef typename Base::difference_type difference_type;
        typedef value_type* pointer;
        typedef value_type& reference;

        VertexIterator();
        VertexIterator(const VertexIterator&);
        VertexIterator(const AdjacencyIterator&);
        VertexIterator(const GraphType&);
        VertexIterator(const GraphType&, const size_type);
        VertexIterator(const GraphType&, const size_type, const size_type);

        // access
        value_type operator*() const;
        value_type operator[](const difference_type) const;

        void coordinate(VertexCoordinate&) const;
    private:
        pointer operator->() const;
    };

    class EdgeIterator
        :          public AdjacencyIterator {
    public:
        typedef GridGraph<DIMENSION, Visitor> GraphType;
        typedef AdjacencyIterator Base;
        typedef const size_type value_type;
        typedef typename Base::difference_type difference_type;
        typedef value_type* pointer;
        typedef value_type& reference;

        EdgeIterator();
        EdgeIterator(const EdgeIterator&);
        EdgeIterator(const AdjacencyIterator&);
        EdgeIterator(const GraphType&);
        EdgeIterator(const GraphType&, const size_type);
        EdgeIterator(const GraphType&, const size_type, const size_type);

        // access
        value_type operator*() const;
        value_type operator[](const difference_type) const;
    private:
        pointer operator->() const;
    };
    // \endcond

    // construction
    GridGraph(const Visitor& = Visitor());
    GridGraph(const VertexCoordinate&, const Visitor& = Visitor());
    GridGraph(const std::initializer_list<std::size_t>, const Visitor& = Visitor());
    void assign(const Visitor& = Visitor());
    void assign(const VertexCoordinate&, const Visitor& = Visitor());

    // iterator access (compatible with Digraph)
    VertexIterator verticesFromVertexBegin(const size_type) const;
    VertexIterator verticesFromVertexEnd(const size_type) const;
    VertexIterator verticesToVertexBegin(const size_type) const;
    VertexIterator verticesToVertexEnd(const size_type) const;
    EdgeIterator edgesFromVertexBegin(const size_type) const;
    EdgeIterator edgesFromVertexEnd(const size_type) const;
    EdgeIterator edgesToVertexBegin(const size_type) const;
    EdgeIterator edgesToVertexEnd(const size_type) const;
    AdjacencyIterator adjacenciesFromVertexBegin(const size_type) const;
    AdjacencyIterator adjacenciesFromVertexEnd(const size_type) const;
    AdjacencyIterator adjacenciesToVertexBegin(const size_type) const;
    AdjacencyIterator adjacenciesToVertexEnd(const size_type) const;

    // access (compatible with Digraph)
    size_type numberOfVertices() const;
    size_type numberOfEdges() const;
    size_type numberOfEdgesFromVertex(const size_type) const;
    size_type numberOfEdgesToVertex(const size_type) const;
    size_type vertexOfEdge(const size_type, const size_type) const;
    size_type edgeFromVertex(const size_type, const size_type) const;
    size_type edgeToVertex(const size_type, const size_type) const;
    size_type vertexFromVertex(const size_type, const size_type) const;
    size_type vertexToVertex(const size_type, const size_type) const;
    AdjacencyType adjacencyFromVertex(const size_type, const size_type) const;
    AdjacencyType adjacencyToVertex(const size_type, const size_type) const;
    std::pair<bool, size_type> findEdge(const size_type, const size_type) const;
    bool multipleEdgesEnabled() const;

    size_type insertEdge(const size_type, const size_type) const;

    size_type shape(const size_type) const;

    size_type vertex(const VertexCoordinate&) const;
    void vertex(const size_type, VertexCoordinate&) const;
    size_type edge(const EdgeCoordinate&) const;
    void edge(const size_type, EdgeCoordinate&) const;

private:
    size_type vertexFromVertex(const VertexCoordinate&, const size_type, size_type&, bool&) const;
    void adjacencyFromVertex(const VertexCoordinate&, const size_type, size_type&, size_type&) const;

    // Member variables
    VertexCoordinate shape_;
    std::array<size_type, DIMENSION> edgeIndexOffsets_;
    std::array<size_type, DIMENSION> vertexIndexOffsets_;
    std::array<VertexCoordinate, DIMENSION> edgeShapes_;
    size_type numberOfVertices_;
    const size_type& numberOfEdges_;
    Visitor visitor_;
};

/// Construct an empty grid graph.
/// \tparam S the type of the indices used. It must be an unsigned type
/// (otherwise a compilation error is caused). Defaults to \c std::size_t .
/// \tparam VISITOR a visitor class to be used. Since this class is immutable
/// appart from resizing, this is a dummy variable meant for compatibility,
/// defaulting to \c andres::graph::IdgeGraphVisitor.
/// \param visitor Visitor to follow changes of integer indices of vertices
/// and edges.
template<unsigned char D, class VISITOR>
inline
GridGraph<D, VISITOR>::GridGraph(
    const Visitor& visitor
)
:   GridGraph(VertexCoordinate({}), visitor) // Chain-call Constructor
{}

/// Construct a grid graph with a specified shape.
///
/// \param shape the shape of the Grid Graph as an std::array.
/// \param visitor Visitor to follow changes of integer indices of vertices
/// and edges.
template<unsigned char D, class VISITOR>
inline
GridGraph<D, VISITOR>::GridGraph(
    const VertexCoordinate& shape,
    const Visitor& visitor
)
    :   numberOfEdges_ (edgeIndexOffsets_[DIMENSION - 1]) {
    assign(shape, visitor);
}

/// Construct a grid graph with a specified shape.
///
/// \param shape the shape of the Grid Graph as an std::array.
/// \param visitor Visitor to follow changes of integer indices of vertices
/// and edges.
template<unsigned char D, class VISITOR>
inline
GridGraph<D, VISITOR>::GridGraph(
    std::initializer_list<std::size_t> shape,
    const Visitor& visitor
)
    :   numberOfEdges_ (edgeIndexOffsets_[DIMENSION - 1]) {
    assert(shape.size()==D);
    VertexCoordinate vCoord;
    std::copy(shape.begin(), shape.end(), vCoord.begin());
    assign(vCoord, visitor);
}

/// Clear a grid graph.
///
/// \param visitor Visitor to follow changes of integer indices of vertices
/// and edges.
///
template<unsigned char D, class VISITOR>
inline void
GridGraph<D, VISITOR>::assign(
    const Visitor& visitor
) {
    VertexCoordinate shape;
    std::fill(shape.begin(), shape.end(), 0);
    assign(shape, visitor);
}

/// Clear a grid graph and assign a new shape.
///
/// \param shape the shape of the grid graph
/// \param visitor Visitor to follow changes of integer indices of vertices
/// and edges.
template<unsigned char D, class VISITOR>
inline void
GridGraph<D, VISITOR>::assign(
    const VertexCoordinate& shape,
    const Visitor& visitor
) {
    shape_ = shape;
    visitor_ = visitor;

    // Set vertex offsets for fast vertex indexing.
    {
        size_type cumprod = 1;
        size_type i;
        for(i = 0; i < DIMENSION; ++i) {
            vertexIndexOffsets_[i] = cumprod;
            cumprod *= shape_[i];
        }
        numberOfVertices_ = cumprod;
    }
    {
        size_type edgeIndexOffset = 0;  // First edge is at offset 0
        for(size_type i = 0; i < DIMENSION; ++i) {
            VertexCoordinate& edgeShape = edgeShapes_[i];
            edgeShape = shape_; // the i-th dimension of the edges along the i-th dimension is 1 less
            if(edgeShape[i] > 0) { // If already zero, no need to reduce.
                --edgeShape[i];
            }
            {
                //
                size_type cumprod = edgeShape[0];
                for(size_type j = 1; j < DIMENSION; ++j) {
                    cumprod *= edgeShape[j];
                }
                edgeIndexOffsets_[i] = (edgeIndexOffset += cumprod);
            }
        }
    }
}

/// Get an iterator to the beginning of the sequence of vertices reachable
/// from a given vertex via a single edge.
///
/// \param vertex Integer index of the vertex.
/// \return VertexIterator.
///
/// \sa verticesFromVertexEnd()
///
template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::VertexIterator
GridGraph<D, VISITOR>::verticesFromVertexBegin(
    const size_type vertex
) const {
    return VertexIterator(*this, vertex, 0);
}

/// Get an iterator to the end of the sequence of vertices reachable from
/// a given vertex via a single edge.
///
/// \param vertex Integer index of the vertex.
/// \return VertexIterator.
///
/// \sa verticesFromVertexBegin()
///
template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::VertexIterator
GridGraph<D, VISITOR>::verticesFromVertexEnd(
    const size_type vertex
) const {
    return VertexIterator(*this, vertex, numberOfEdgesFromVertex(vertex));
}

/// Get an iterator to the beginning of the sequence of vertices from which
/// a given vertex is reachable via a single edge.
///
/// \param vertex Integer index of the vertex.
/// \return VertexIterator.
///
/// \sa verticesToVertexEnd()
///
template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::VertexIterator
GridGraph<D, VISITOR>::verticesToVertexBegin(
    const size_type vertex
) const {
    return VertexIterator(*this, vertex, 0);
}

/// Get an iterator to the end of the sequence of vertices from which a
/// given vertex is reachable via a single edge.
///
/// \param vertex Integer index of the vertex.
/// \return VertexIterator.
///
/// \sa verticesToVertexBegin()
///
template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::VertexIterator
GridGraph<D, VISITOR>::verticesToVertexEnd(
    const size_type vertex
) const {
    return VertexIterator(*this, vertex, numberOfEdgesFromVertex(vertex));
}

/// Get an iterator to the beginning of the sequence of edges that originate
/// from a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return EdgeIterator.
///
/// \sa edgesFromVertexEnd()
///
template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::EdgeIterator
GridGraph<D, VISITOR>::edgesFromVertexBegin(
    const size_type vertex
) const {
    return EdgeIterator(*this, vertex, 0);
}

/// Get an iterator to the end of the sequence of edges that originate from
/// a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return EdgeIterator.
///
/// \sa edgesFromVertexBegin()
///
template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::EdgeIterator
GridGraph<D, VISITOR>::edgesFromVertexEnd(
    const size_type vertex
) const {
    return EdgeIterator(*this, vertex, numberOfEdgesFromVertex(vertex));
}

/// Get an iterator to the beginning of the sequence of edges that are
/// incident to a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return EdgeIterator.
///
/// \sa edgesToVertexEnd()
///
template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::EdgeIterator
GridGraph<D, VISITOR>::edgesToVertexBegin(
    const size_type vertex
) const {
    return EdgeIterator(*this, vertex, 0);
}

/// Get an iterator to the end of the sequence of edges that are incident
/// to a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return EdgeIterator.
///
/// \sa edgesToVertexBegin()
///
template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::EdgeIterator
GridGraph<D, VISITOR>::edgesToVertexEnd(
    const size_type vertex
) const {
    return EdgeIterator(*this, vertex, numberOfEdgesFromVertex(vertex));
}

/// Get an iterator to the beginning of the sequence of adjacencies that
/// originate from a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return AdjacencyIterator.
///
/// \sa adjacenciesFromVertexEnd()
///
template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::AdjacencyIterator
GridGraph<D, VISITOR>::adjacenciesFromVertexBegin(
    const size_type vertex
) const {
    return AdjacencyIterator(*this, vertex, 0);
}

/// Get an iterator to the end of the sequence of adjacencies that originate
/// from a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return AdjacencyIterator.
///
/// \sa adjacenciesFromVertexBegin()
///
template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::AdjacencyIterator
GridGraph<D, VISITOR>::adjacenciesFromVertexEnd(
    const size_type vertex
) const {
    return AdjacencyIterator(*this, vertex, numberOfEdgesFromVertex(vertex));
}

/// Get an iterator to the beginning of the sequence of adjacencies incident
/// to a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return AdjacencyIterator.
///
/// \sa adjacenciesToVertexEnd()
///
template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::AdjacencyIterator
GridGraph<D, VISITOR>::adjacenciesToVertexBegin(
    const size_type vertex
) const {
    return AdjacencyIterator(*this, vertex, 0);
}

/// Get an iterator to the end of the sequence of adjacencies incident to
/// a given vertex.
///
/// \param vertex Integer index of the vertex.
/// \return AdjacencyIterator.
///
/// \sa adjacenciesToVertexBegin()
///
template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::AdjacencyIterator
GridGraph<D, VISITOR>::adjacenciesToVertexEnd(
    const size_type vertex
) const {
    return AdjacencyIterator(*this, vertex, numberOfEdgesFromVertex(vertex));
}

/// Get the number of vertices.
///
template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::size_type
GridGraph<D, VISITOR>::numberOfVertices() const {
    return numberOfVertices_;
}

/// Get the number of edges.
///
template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::size_type
GridGraph<D, VISITOR>::numberOfEdges() const {
    return numberOfEdges_;
}

/// Get the number of edges that originate from a given vertex.
///
/// \param vertex Integer index of a vertex.
///
/// \sa edgeFromVertex()
///
template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::size_type
GridGraph<D, VISITOR>::numberOfEdgesFromVertex(
    const size_type vertex
) const {
    assert(vertex < numberOfVertices());
    VertexCoordinate edgeCoordinate;
    this->vertex(vertex, edgeCoordinate);
    size_type numEdgesFromVertex = 0;
    for(size_type i = 0; i < DIMENSION; ++i) {
        const size_type& coordinate = edgeCoordinate[i];
        if(coordinate > 0) {
            ++numEdgesFromVertex;
        }
        if(coordinate < (shape_[i] - 1)) {
            ++numEdgesFromVertex;
        }
    }
    return numEdgesFromVertex;
}

/// Get the number of  edges that are incident to a given vertex.
///
/// \param vertex Integer index of a vertex.
///
/// \sa edgeToVertex()
///
template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::size_type
GridGraph<D, VISITOR>::numberOfEdgesToVertex(
    const size_type vertex
) const {
    return numberOfEdgesFromVertex(vertex);
}

/// Get the integer index of a vertex of an edge.
///
/// \param edge Integer index of an edge.
/// \param j Number of the vertex in the edge; either 0 or 1.
///
template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::size_type
GridGraph<D, VISITOR>::vertexOfEdge(
    const size_type edge,
    const size_type j
) const {
    assert(edge < numberOfEdges());
    assert(j < 2);
    EdgeCoordinate edgeCoordinate;
    this->edge(edge, edgeCoordinate);
    size_type pivotIndex = vertex(edgeCoordinate.pivotCoordinate_);
    if(j == 0) {
        return pivotIndex;
    }
    else {
        return pivotIndex + vertexIndexOffsets_[edgeCoordinate.dimension_];
    }
}

/// Get the integer index of an edge that originates from a given vertex.
///
/// \param vertex Integer index of a vertex.
/// \param j Number of the edge; between 0 and numberOfEdgesFromVertex(vertex)
/// - 1.
///
/// \sa numberOfEdgesFromVertex()exFi
///
template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::size_type
GridGraph<D, VISITOR>::edgeFromVertex(
    const size_type vertex,
    const size_type j
) const {
    assert(j < numberOfEdgesFromVertex(vertex));
    VertexCoordinate vertexCoordinate;
    this->vertex(vertex, vertexCoordinate);
    size_type direction;
    bool isSmaller;
    vertexFromVertex(vertexCoordinate, j, direction, isSmaller);
    if(isSmaller) {
        VertexCoordinate pivotVertexCoordinate = vertexCoordinate;
        --pivotVertexCoordinate[direction];
        return edge(EdgeCoordinate(pivotVertexCoordinate, direction));
    }
    else {
        return edge(EdgeCoordinate(vertexCoordinate, direction));
    }
}

/// Get the integer index of an edge that is incident to a given vertex.
///
/// \param vertex Integer index of a vertex.
/// \param j Number of the edge; between 0 and numberOfEdgesFromVertex(vertex)
/// - 1.
///
/// \sa numberOfEdgesToVertex()
///
template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::size_type
GridGraph<D, VISITOR>::edgeToVertex(
    const size_type vertex,
    const size_type j
) const {
    assert(j < numberOfEdgesToVertex(vertex));
    return edgeFromVertex(vertex, j);
}

/// Get the integer index of a vertex from which a given vertex is reachable
/// via a single edge.
///
/// \param vertex Integer index of a vertex.
/// \param j Number of the vertex; between 0 and
/// numberOfEdgesFromVertex(vertex) - 1.
///
/// \sa numberOfEdgesFromVertex()
///
template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::size_type
GridGraph<D, VISITOR>::vertexFromVertex(
    const size_type vertex,
    const size_type j
) const {
    assert(j < numberOfEdgesToVertex(vertex));
    VertexCoordinate vertexCoordinate;
    this->vertex(vertex, vertexCoordinate);
    size_type direction;
    bool isSmaller;
    return vertexFromVertex(vertexCoordinate, j, direction, isSmaller);
}


/// Get the integer index of a vertex to which a given vertex has a single
/// edge.
///
/// \param vertex Integer index of a vertex.
/// \param j Number of the vertex; between 0 and
/// numberOfEdgesFromVertex(vertex) - 1.
///
/// \sa numberOfEdgesFromVertex()
///
template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::size_type
GridGraph<D, VISITOR>::vertexToVertex(
    const size_type vertex,
    const size_type j
) const {
    assert(j < numberOfEdgesToVertex(vertex));
    return vertexFromVertex(vertex, j);
}

/// Get the j-th adjacency from a vertex.
///
/// \param vertex Vertex.
/// \param j Number of the adjacency.
template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::AdjacencyType
GridGraph<D, VISITOR>::adjacencyFromVertex(
    const size_type vertex,
    const size_type j
) const {
    assert(j < numberOfEdgesToVertex(vertex));
    size_type direction;
    size_type adjacentEdgeIndex;
    size_type adjacentVertexIndex;

    {
        VertexCoordinate vertexCoordinate;
        this->vertex(vertex, vertexCoordinate);

        adjacencyFromVertex(vertexCoordinate, j, adjacentVertexIndex, adjacentEdgeIndex);
    }
    return AdjacencyType(adjacentVertexIndex, adjacentEdgeIndex);
}


/// Get the j-th adjacency to a vertex.
///
/// \param vertex Vertex.
/// \param j Number of the adjacency.
///
template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::AdjacencyType
GridGraph<D, VISITOR>::adjacencyToVertex(
    const size_type vertex,
    const size_type j
) const {
    return adjacencyFromVertex(vertex, j);
}

/// Search for an edge (in constant time).
///
/// \param vertex0 first vertex of the edge.
/// \param vertex1 second vertex of the edge.
/// \retval pair an \c std::pair. If an edge from \b vertex0 to \b vertex1
/// exists, \c pair.first is \c true
/// and \c pair.second is the index of such an edge. If no edge from \b
/// vertex0 to \b vertex1 exists, \c pair.first is \c false
/// and \c pair.second is undefined.
template<unsigned char D, class VISITOR>
inline std::pair<bool, typename GridGraph<D, VISITOR>::size_type>
GridGraph<D, VISITOR>::findEdge(
    const size_type vertex0,
    const size_type vertex1
) const {
    assert(vertex0 < numberOfVertices());
    assert(vertex1 < numberOfVertices());
    if(vertex0 < vertex1) {
        size_type offset = vertex1 - vertex0;
        size_type i;
        for(i = 0; i < DIMENSION - 1; ++i) {
            if(offset == vertexIndexOffsets_[i]) { // edge found unless offset runs across dimension
                VertexCoordinate vertexCoordinate0;
                vertex(vertex0, vertexCoordinate0);
                if(vertexCoordinate0[i] != shape_[i] - 1) { // Check for boundary case
                    const EdgeCoordinate edgeCoordinate(vertexCoordinate0, i, false);
                    const size_type edgeIndex = edge(edgeCoordinate);
                    return std::make_pair(true, edgeIndex);
                }
            }
        }
        if(offset == vertexIndexOffsets_[i]) { // edge found
            VertexCoordinate vertexCoordinate0;
            vertex(vertex0, vertexCoordinate0);
            const EdgeCoordinate edgeCoordinate(vertexCoordinate0, i, false);
            const size_type edgeIndex = edge(edgeCoordinate);
            return std::make_pair(true, edgeIndex);
        }
    }
    else {   // On expectation faster to ignore the equal case. (Assuming uniform input distribution)
        size_type offset = vertex0 - vertex1;
        size_type i;
        for(i = 0; i < DIMENSION - 1; ++i) {
            if(offset == vertexIndexOffsets_[i]) { // edge found unless offset runs across dimensions
                VertexCoordinate vertexCoordinate1;
                vertex(vertex1, vertexCoordinate1);
                if(vertexCoordinate1[i] != shape_[i] - 1) { // Check for boundary case
                    const EdgeCoordinate edgeCoordinate(vertexCoordinate1, i, false);
                    const size_type edgeIndex = edge(edgeCoordinate);
                    return std::make_pair(true, edgeIndex);
                }
            }
        }
        if(offset == vertexIndexOffsets_[i]) { // edge found unless offset runs across dimensions
            VertexCoordinate vertexCoordinate1;
            vertex(vertex1, vertexCoordinate1);
            const EdgeCoordinate edgeCoordinate(vertexCoordinate1, i, false);
            const size_type edgeIndex = edge(edgeCoordinate);
            return std::make_pair(true, edgeIndex);
        }
    }
    return std::make_pair(false, 0);
}

/// Indicate if multiple edges are enabled.
///
/// \return false
///
template<unsigned char D, class VISITOR>
inline bool
GridGraph<D, VISITOR>::multipleEdgesEnabled() const {
    return false;
}

/// Returns the edge between the specified vertices.
///
/// \param vertexIndex0 Integer index of the first vertex in the edge.
/// \param vertexIndex1 Integer index of the second vertex in the edge.
/// \return Integer index of the newly inserted edge.
/// \throw runtime_error If the edge does not exist.
template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::size_type
GridGraph<D, VISITOR>::insertEdge(
    const size_type vertexIndex0,
    const size_type vertexIndex1
) const {
    assert(vertexIndex0 < numberOfVertices());
    assert(vertexIndex1 < numberOfVertices());

    std::pair<bool, std::size_t> p = findEdge(vertexIndex0, vertexIndex1);
    if(p.first == true) {
        return p.second;
    } else {
        throw std::runtime_error("Edge not found.");
    }
}

/// Retrieve the specified vertex of the graph.
/// \param[in] vertexIndex the integer index of the requested vertex
/// \param[out] vertexCoordinate The coordinates of the vertex.
/// \warning For the sake of performance this function does not validate
/// its inputs.
template<unsigned char D, class VISITOR>
void
GridGraph<D, VISITOR>::vertex(
    size_type vertexIndex,
    VertexCoordinate& vertexCoordinate
) const {
    assert(vertexIndex < numberOfVertices_);
    size_type i;
    for(i = 0; i < DIMENSION - 1; ++i) {
        vertexCoordinate[i] = vertexIndex % shape_[i];
        vertexIndex = vertexIndex / shape_[i];
    }
    vertexCoordinate[i] = vertexIndex;
}

/// Get the size of a specific dimension of the grid graph.
/// \param dimension the index of the dimension index to retrieve.
/// \return the size of the specified dimension.
/// \see mapIndexToCoordinate mapVertexCoordinateToIndex
template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::size_type
GridGraph<D, VISITOR>::shape(
    const size_type dimension
) const {
    return shape_[dimension];
}

/// Retrieve the specified vertex of the graph.
/// \param vertexCoordinate the integer index of the requested vertex
/// \return The integer index of the specified vertex.
/// \warning For the sake of performance this function does not validate
/// its inputs.
template<unsigned char D, class VISITOR>
typename GridGraph<D, VISITOR>::size_type
GridGraph<D, VISITOR>::vertex(
    const VertexCoordinate& vertexCoordinate
) const {
    size_type index = vertexCoordinate[DIMENSION - 1];
    for(size_type i = DIMENSION - 1; i > 0; --i) {
        index = index * shape_[i - 1] + vertexCoordinate[i - 1];
    }
    return index;
}

/// Retrieve the specified edge of the graph.
/// \param edgeCoordinate the coordinates of the minimum endpoint (\e pivot)
/// and the direction of the requested edge.
/// \return The integer index of the specified edge.
/// \warning For the sake of performance this function does not validate
/// its inputs.
/// \sa hasEdge()
template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::size_type
GridGraph<D, VISITOR>::edge(
    const EdgeCoordinate& edgeCoordinate
) const {
    assert(edgeCoordinate.dimension_ < DIMENSION);
    const size_type& dimension = edgeCoordinate.dimension_;
    const VertexCoordinate& pivotCoordinate = edgeCoordinate.pivotCoordinate_;
    const VertexCoordinate& edgeShape = edgeShapes_[dimension];

    // size_type index = mapCoordinateToIndex(edgeCoordinate,edgeShape);
    size_type index = pivotCoordinate[DIMENSION - 1];
    for(size_type i = DIMENSION - 1; i > 0; --i) {
        index = index * edgeShape[i - 1] + pivotCoordinate[i - 1];
    }
    if(dimension > 0) {
        const size_type& offset = edgeIndexOffsets_[dimension - 1];
        index += offset;
    }
    return index;
}

/// Retrieve the specified edge of the graph.
/// \param[in] edgeIndex the integer index of the requested edge.
/// \param[out] edgeCoordinate a GridGraph::EdgeCoordinate instance. \c
/// edgeCoordinate.pivot is the integer index
///  of the minimum of the two edge endpoints; \p edgeCoordinate.direction
///  is the dimension along which
///  the edge is drawn, with an assumed positive isSmaller.
/// \warning For the sake of performance this function does not validate
/// its inputs.
/// \sa GridGraph::bool
template<unsigned char D, class VISITOR>
inline void
GridGraph<D, VISITOR>::edge(
    size_type edgeIndex,
    EdgeCoordinate& edgeCoordinate
) const {
    // WARNING: If the assertion does not hold, the code will scan until an unlikely condition.
    assert(edgeIndex < numberOfEdges());
    size_type& direction = edgeCoordinate.dimension_;
    // Find the direction as the last edge offset:
    for(direction = 0; edgeIndex >= edgeIndexOffsets_[direction]; ++direction);
    if(direction > 0) { // Not needed, but saves one memory lookup.
        const size_type& offset = edgeIndexOffsets_[direction - 1];
        edgeIndex -= offset;
    }
    const VertexCoordinate& edgeShape = edgeShapes_[direction];
    // mapEdgeIndexToCoordinate
    {
        VertexCoordinate& pivotCoordinate = edgeCoordinate.pivotCoordinate_;
        size_type i;
        for(i = 0; i < DIMENSION - 1; ++i) {
            pivotCoordinate[i] = edgeIndex % edgeShape[i];
            edgeIndex = edgeIndex / edgeShape[i];
        }
        pivotCoordinate[i] = edgeIndex;
    }
}


/// \internal
/// \brief Retrieve the direction and isSmaller of the <b>j</b>-th edge of
/// a vertex
/// \param vertexCoordinate the integer index of the vertex from which the
/// edges are originating
/// \param j the index within the set of adacent edges of the specified vertex
/// \param[out] direction the direction of the specified edge
/// \param[out] isSmaller the isSmaller of the specified edge
/// \retval found If the edge is found, a value of \c true is returned.
/// \remark This function attempts a fast imlementation that tries not
/// to waste any comparison operations and is meant for use as a building
/// block of the class.
/// \warning For the sake of performance this function does not validate
/// its inputs.
/// \sa directionOfEdgeFromVertex, edgeFromVertex
template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::size_type
GridGraph<D, VISITOR>::vertexFromVertex(
    const VertexCoordinate& vertexCoordinate,
    const size_type j,
    size_type& direction,
    bool& isSmaller
) const {
    assert(vertex(vertexCoordinate) < numberOfVertices());
    VertexCoordinate modifieableVertexCoordinate = vertexCoordinate;
    size_type cur_j = 0;
    for(size_type i = DIMENSION; i > 0; --i) {
        if(vertexCoordinate[i - 1] > 0) { // edge exists
            if(cur_j == j) {
                direction = i - 1;
                isSmaller = true;
                --modifieableVertexCoordinate[i - 1];
                const size_type vertexIndex = vertex(modifieableVertexCoordinate);
                return vertexIndex;
            }
            ++cur_j;
        }
    }
    for(size_type i = 0; i < DIMENSION; ++i) {
        if(vertexCoordinate[i] < shape_[i] - 1) { // edge exists
            if(cur_j == j) {
                direction = i;
                isSmaller = false;
                ++modifieableVertexCoordinate[i];
                const size_type vertexIndex = vertex(modifieableVertexCoordinate);
                return vertexIndex;
            }
            ++cur_j;
        }
    }
    throw std::out_of_range("vertex neighbor index out of range.");
}

template<unsigned char D, class VISITOR>
inline void
GridGraph<D, VISITOR>::adjacencyFromVertex(
    const VertexCoordinate& vertexCoordinate,
    const size_type j,
    size_type& adjacentVertexIndex,
    size_type& adjacentEdgeIndex
) const {
    assert(vertex(vertexCoordinate) < numberOfVertices());
    size_type direction;
    bool isSmaller;
    adjacentVertexIndex = vertexFromVertex(vertexCoordinate, j, direction, isSmaller);
    if(isSmaller) {
        VertexCoordinate pivotVertexCoordinate = vertexCoordinate;
        --pivotVertexCoordinate[direction];
        adjacentEdgeIndex = edge(EdgeCoordinate(pivotVertexCoordinate, direction));
    }
    else {
        adjacentEdgeIndex = edge(EdgeCoordinate(vertexCoordinate, direction));
    }

}

/// Initialize an edge coordinate.
/// \param pivotCoordinate coordinate of the reference vertex.
/// \param dimension dimension along which the edge is drawn.
/// \param isSmaller relative position of specified endpoint.\n
/// A value of \c true results in an object corresponding to the edge
/// of which the smallest endpoint has a GridGraph::VertexCoordinate
/// equal to \b pivotCoordinate.\n
/// Similarly, a value of \c false results in an object corresponding
/// to the edge of which the largest endpoint has a
/// GridGraph::VertexCoordinate equal to \b pivotCoordinate.
template<unsigned char D, class VISITOR>
inline GridGraph<D, VISITOR>::EdgeCoordinate::EdgeCoordinate(
    const VertexCoordinate& pivotCoordinate,
    const size_type dimension,
    const bool isSmaller
)
    :   pivotCoordinate_(pivotCoordinate),
        dimension_(dimension) {
    if(isSmaller) {
        assert(pivotCoordinate_[dimension] > 0);
        --pivotCoordinate_[dimension];
    }
}

/// Default non-initializing constructor
/// \warning this constructor will \b NOT set the pivotCoordinate_ membr
/// variable to zero.
template<unsigned char D, class VISITOR>
inline GridGraph<D, VISITOR>::EdgeCoordinate::EdgeCoordinate() {
}


// \cond SUPPRESS_DOXYGEN

// implementation of AdjacencyIterator

template<unsigned char D, class VISITOR>
inline
GridGraph<D, VISITOR>::AdjacencyIterator::AdjacencyIterator()
    :   vertex_(0),
        adjacencyIndex_(0),
        adjacency_()
{}
template<unsigned char D, class VISITOR>
inline
GridGraph<D, VISITOR>::AdjacencyIterator::AdjacencyIterator(
    const GraphType& graph
)
    :   graph_(&graph),
        vertex_(0),
        adjacencyIndex_(0),
        adjacency_()
{}

template<unsigned char D, class VISITOR>
inline
GridGraph<D, VISITOR>::AdjacencyIterator::AdjacencyIterator(
    const GraphType& graph,
    const size_type vertex
)
    :   graph_(&graph),
        vertex_(vertex),
        adjacencyIndex_(0),
        adjacency_() {
    assert(vertex < graph.numberOfVertices());
}

template<unsigned char D, class VISITOR>
inline
GridGraph<D, VISITOR>::AdjacencyIterator::AdjacencyIterator(
    const GraphType& graph,
    const size_type vertex,
    const size_type adjacencyIndex
)
    :   graph_(&graph),
        vertex_(vertex),
        adjacencyIndex_(adjacencyIndex),
        adjacency_() {
    assert(vertex < graph.numberOfVertices());
    assert(adjacencyIndex <= graph.numberOfVertices());
}

template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::AdjacencyIterator&
GridGraph<D, VISITOR>::AdjacencyIterator::operator+=(
    const difference_type d
) {
    adjacencyIndex_ += d;
    return *this;
}

template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::AdjacencyIterator&
GridGraph<D, VISITOR>::AdjacencyIterator::operator-=(
    const difference_type d
) {
    adjacencyIndex_ -= d;
    return *this;
}

template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::AdjacencyIterator&
GridGraph<D, VISITOR>::AdjacencyIterator::operator++() {
    ++adjacencyIndex_;
    return *this;
}

template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::AdjacencyIterator&
GridGraph<D, VISITOR>::AdjacencyIterator::operator--() {
    --adjacencyIndex_;
    return *this;
}

template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::AdjacencyIterator
GridGraph<D, VISITOR>::AdjacencyIterator::operator++(int) {
    AdjacencyIterator copy = *this;
    ++adjacencyIndex_;
    return copy;
}

template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::AdjacencyIterator
GridGraph<D, VISITOR>::AdjacencyIterator::operator--(int) {
    AdjacencyIterator copy = *this;
    --adjacencyIndex_;
    return copy;
}

template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::AdjacencyIterator
GridGraph<D, VISITOR>::AdjacencyIterator::operator+(
    const difference_type d
) const {
    AdjacencyIterator copy = *this;
    copy += d;
    return copy;
}

template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::AdjacencyIterator
GridGraph<D, VISITOR>::AdjacencyIterator::operator-(
    const difference_type d
) const {
    AdjacencyIterator copy = *this;
    copy -= d;
    return copy;
}

template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::AdjacencyIterator::difference_type
GridGraph<D, VISITOR>::AdjacencyIterator::operator-(
    const AdjacencyIterator& adjacencyIterator
) const {
    return adjacencyIndex_ - adjacencyIterator.adjacencyIndex_;
}

template<unsigned char D, class VISITOR>
inline bool
GridGraph<D, VISITOR>::AdjacencyIterator::operator==(
    const AdjacencyIterator& other
) const {
    return adjacencyIndex_ == other.adjacencyIndex_
           && vertex_ == other.vertex_
           && graph_ == other.graph_;
}

template<unsigned char D, class VISITOR>
inline bool
GridGraph<D, VISITOR>::AdjacencyIterator::operator!=(
    const AdjacencyIterator& other
) const {
    return adjacencyIndex_ != other.adjacencyIndex_
           || vertex_ != other.vertex_
           || graph_ != other.graph_;
}

template<unsigned char D, class VISITOR>
inline bool
GridGraph<D, VISITOR>::AdjacencyIterator::operator<(
    const AdjacencyIterator& other
) const {
    return adjacencyIndex_ < other.adjacencyIndex_
           && vertex_ == other.vertex_
           && graph_ == other.graph_;
}

template<unsigned char D, class VISITOR>
inline bool
GridGraph<D, VISITOR>::AdjacencyIterator::operator<=(
    const AdjacencyIterator& other
) const {
    return adjacencyIndex_ <= other.adjacencyIndex_
           && vertex_ == other.vertex_
           && graph_ == other.graph_;
}

template<unsigned char D, class VISITOR>
inline bool
GridGraph<D, VISITOR>::AdjacencyIterator::operator>(
    const AdjacencyIterator& other
) const {
    return adjacencyIndex_ > other.adjacencyIndex_
           && vertex_ == other.vertex_
           && graph_ == other.graph_;
}

template<unsigned char D, class VISITOR>
inline bool
GridGraph<D, VISITOR>::AdjacencyIterator::operator>=(
    const AdjacencyIterator& other
) const {
    return adjacencyIndex_ >= other.adjacencyIndex_
           && vertex_ == other.vertex_
           && graph_ == other.graph_;
}

// Since the GridGraph has no backing storage for each of the Adjacency
// objects handled in the class, the AdjacencyIterator is limited to only
// one adjacency object per iterator instance.
// This should be sufficient for most normal usage scenarios, and is supported
// by the very lightweight copying of the AdjacencyType object itself.
// Note, however, that if you store a reference to the pointed object,
// then advance the iterator and subsequently dereference it again,
// the new reference will refer to the very same same adjacency object,
// unique for the AdjacencyIterator instance.
// This operation will therefore silently update all previous references to
// the adjacency object of the iterator.
template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::AdjacencyIterator::reference
GridGraph<D, VISITOR>::AdjacencyIterator::operator*() {
    adjacency_ = graph_->adjacencyFromVertex(vertex_, adjacencyIndex_);
    return adjacency_;
}

template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::AdjacencyIterator::pointer
GridGraph<D, VISITOR>::AdjacencyIterator::operator->() {
    adjacency_ = graph_->adjacencyFromVertex(vertex_, adjacencyIndex_);
    return &adjacency_;
}

template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::AdjacencyIterator::reference
GridGraph<D, VISITOR>::AdjacencyIterator::operator[](
    const difference_type j
) {
    adjacency_ = graph_->adjacencyFromVertex(vertex_, adjacencyIndex_ + j);
    return adjacency_;
}

// implementation of VertexIterator

template<unsigned char D, class VISITOR>
inline
GridGraph<D, VISITOR>::VertexIterator::VertexIterator()
    :   Base()
{}

template<unsigned char D, class VISITOR>
inline
GridGraph<D, VISITOR>::VertexIterator::VertexIterator(
    const GraphType& graph
)
    :   Base(graph)
{}

template<unsigned char D, class VISITOR>
inline
GridGraph<D, VISITOR>::VertexIterator::VertexIterator(
    const GraphType& graph,
    const size_type vertex
)
    :   Base(graph, vertex)
{}

template<unsigned char D, class VISITOR>
inline
GridGraph<D, VISITOR>::VertexIterator::VertexIterator(
    const GraphType& graph,
    const size_type vertex,
    const size_type adjacencyIndex
)
    :   Base(graph, vertex, adjacencyIndex)
{}

template<unsigned char D, class VISITOR>
inline
GridGraph<D, VISITOR>::VertexIterator::VertexIterator(
    const VertexIterator& it
)
    :   Base(*(it.graph_), it.vertex_, it.adjacencyIndex_)
{}

template<unsigned char D, class VISITOR>
inline
GridGraph<D, VISITOR>::VertexIterator::VertexIterator(
    const AdjacencyIterator& it
)
    :   Base(it)
{}

template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::VertexIterator::value_type
GridGraph<D, VISITOR>::VertexIterator::operator*() const {
    return Base::graph_->vertexFromVertex(Base::vertex_, Base::adjacencyIndex_);
}

template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::VertexIterator::value_type
GridGraph<D, VISITOR>::VertexIterator::operator[](
    const difference_type j
) const {
    return Base::graph_->vertexFromVertex(Base::vertex_, Base::adjacencyIndex_ + j);
}

template<unsigned char D, class VISITOR>
inline void
GridGraph<D, VISITOR>::VertexIterator::coordinate(
    VertexCoordinate& vertexCoordinate
) const {
    const size_type opposite = Base::graph_->vertexFromVertex(Base::vertex_, Base::adjacencyIndex_);
    Base::graph_->vertex(opposite, vertexCoordinate);
}

// implementation of EdgeIterator

template<unsigned char D, class VISITOR>
inline
GridGraph<D, VISITOR>::EdgeIterator::EdgeIterator()
    :   Base()
{}

template<unsigned char D, class VISITOR>
inline
GridGraph<D, VISITOR>::EdgeIterator::EdgeIterator(
    const GraphType& graph
)
    :   Base(graph)
{}

template<unsigned char D, class VISITOR>
inline
GridGraph<D, VISITOR>::EdgeIterator::EdgeIterator(
    const GraphType& graph,
    const size_type vertex
)
    :   Base(graph, vertex)
{}

template<unsigned char D, class VISITOR>
inline
GridGraph<D, VISITOR>::EdgeIterator::EdgeIterator(
    const GraphType& graph,
    const size_type vertex,
    const size_type adjacencyIndex
)
    :   Base(graph, vertex, adjacencyIndex)
{}

template<unsigned char D, class VISITOR>
inline
GridGraph<D, VISITOR>::EdgeIterator::EdgeIterator(
    const EdgeIterator& it
)
    :   Base(*(it.graph_), it.vertex_, it.adjacencyIndex_)
{}

template<unsigned char D, class VISITOR>
inline
GridGraph<D, VISITOR>::EdgeIterator::EdgeIterator(
    const AdjacencyIterator& it
)
    :   Base(it)
{}

template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::EdgeIterator::value_type
GridGraph<D, VISITOR>::EdgeIterator::operator*() const {
    return Base::graph_->edgeFromVertex(Base::vertex_, Base::adjacencyIndex_);
}

template<unsigned char D, class VISITOR>
inline typename GridGraph<D, VISITOR>::EdgeIterator::value_type
GridGraph<D, VISITOR>::EdgeIterator::operator[](
    const difference_type j
) const {
    return Base::graph_->edgeFromVertex(Base::vertex_, Base::adjacencyIndex_ + j);
}

// \endcond


} // namespace graph
} // namespace andres

#endif // #ifndef ANDRES_GRAPH_GRID_GRAPH_HXX


