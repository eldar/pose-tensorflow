#pragma once
#ifndef ANDRES_GRAPH_ADJACENCY_HXX
#define ANDRES_GRAPH_ADJACENCY_HXX

namespace andres {
namespace graph {

/// The adjacency of a vertex consists of a vertex and a connecting edge.
template<class T = std::size_t>
class Adjacency {
public:
    typedef T Value;

    Adjacency(const Value = T(), const Value = T());
    Value vertex() const;
    Value& vertex();
    Value edge() const;
    Value& edge();
    bool operator<(const Adjacency<Value>&) const;
    bool operator<=(const Adjacency<Value>&) const;
    bool operator>(const Adjacency<Value>&) const;
    bool operator>=(const Adjacency<Value>&) const;
    bool operator==(const Adjacency<Value>&) const;
    bool operator!=(const Adjacency<Value>&) const;

private:
    Value vertex_;
    Value edge_;
};

/// Construct an adjacency.
///
/// \param vertex Vertex.
/// \param edge Edge.
///
template<class T>
inline
Adjacency<T>::Adjacency(
    const Value vertex,
    const Value edge
)
:   vertex_(vertex),
    edge_(edge)
{}

/// Access the vertex.
///
template<class T>
inline typename Adjacency<T>::Value
Adjacency<T>::vertex() const {
    return vertex_;
}

/// Access the vertex.
///
template<class T>
inline typename Adjacency<T>::Value&
Adjacency<T>::vertex() {
    return vertex_;
}

/// Access the edge.
///
template<class T>
inline typename Adjacency<T>::Value
Adjacency<T>::edge() const {
    return edge_;
}

/// Access the edge.
///
template<class T>
inline typename Adjacency<T>::Value&
Adjacency<T>::edge() {
    return edge_;
}

/// Adjacencies are ordered first wrt the vertex, then wrt the edge.
///
template<class T>
inline bool
Adjacency<T>::operator<(
    const Adjacency<T>& in
) const {
    if(vertex_ < in.vertex_) {
        return true;
    }
    else if(vertex_ == in.vertex_) {
        return edge_ < in.edge_;
    }
    else {
        return false;
    }
}

/// Adjacencies are ordered first wrt the vertex, then wrt the edge.
///
template<class T>
inline bool
Adjacency<T>::operator<=(
    const Adjacency<T>& in
) const {
    if(vertex_ < in.vertex_) {
        return true;
    }
    else if(vertex_ == in.vertex_) {
        return edge_ <= in.edge_;
    }
    else {
        return false;
    }
}

/// Adjacencies are ordered first wrt the vertex, then wrt the edge.
///
template<class T>
inline bool
Adjacency<T>::operator>(
    const Adjacency<T>& in
) const {
    if(vertex_ > in.vertex_) {
        return true;
    }
    else if(vertex_ == in.vertex_) {
        return edge_ > in.edge_;
    }
    else {
        return false;
    }
}

/// Adjacencies are ordered first wrt the vertex, then wrt the edge.
///
template<class T>
inline bool
Adjacency<T>::operator>=(
    const Adjacency<T>& in
) const {
    if(vertex_ > in.vertex_) {
        return true;
    }
    else if(vertex_ == in.vertex_) {
        return edge_ >= in.edge_;
    }
    else {
        return false;
    }
}

/// Adjacencies are equal if both the vertex and the edge are equal.
///
template<class T>
inline bool
Adjacency<T>::operator==(
    const Adjacency<T>& in
) const {
    return vertex_ == in.vertex_ && edge_ == in.edge_;
}

/// Adjacencies are unequal if either the vertex or the edge are unqual.
///
template<class T>
inline bool
Adjacency<T>::operator!=(
    const Adjacency<T>& in
) const {
    return !(*this == in);
}

} // namespace graph
} // namespace andres

#endif // #ifndef ANDRES_GRAPH_ADJACENCY_HXX
