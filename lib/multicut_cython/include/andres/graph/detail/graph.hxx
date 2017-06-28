#pragma once
#ifndef ANDRES_GRAPH_DETAIL_GRAPH_HXX
#define ANDRES_GRAPH_DETAIL_GRAPH_HXX

#include <iterator>
#include <andres/random-access-set.hxx>

namespace andres {
namespace graph {
// \cond SUPPRESS_DOXYGEN
namespace detail {

template<bool DIRECTED, class T = std::size_t>
class Edge {
public:
    typedef T Value;

    Edge(const Value, const Value);
    Value operator[](const std::size_t) const;
    Value& operator[](const std::size_t);

private:
    Value vertexIndices_[2];
};

typedef RandomAccessSet<Adjacency<> > Adjacencies;

template<bool T>
class IteratorHelper
:   public Adjacencies::const_iterator
{
private:
    typedef typename Adjacencies::const_iterator Base;

public:
    typedef typename Base::iterator_category iterator_category;
    typedef typename Base::difference_type difference_type;
    typedef const std::size_t value_type;
    typedef value_type* pointer;
    typedef value_type& reference;

    // construction and assignment
    IteratorHelper();
    IteratorHelper(const Base&);
    IteratorHelper(const IteratorHelper<T>&);
    IteratorHelper operator=(const Base&);
    IteratorHelper operator=(const IteratorHelper<T>&);

    // increment and decrement
    IteratorHelper<T>& operator+=(const difference_type);
    IteratorHelper<T>& operator-=(const difference_type);
    IteratorHelper<T>& operator++(); // prefix
    IteratorHelper<T>& operator--(); // prefix
    IteratorHelper<T> operator++(int); // postfix
    IteratorHelper<T> operator--(int); // postfix
    IteratorHelper<T> operator+(const difference_type) const;
    IteratorHelper<T> operator-(const difference_type) const;
    #ifdef _MSC_VER
    difference_type operator-(const IteratorHelper<T>&) const;
    #endif

    // access
    value_type operator*() const;
    value_type operator[](const std::size_t j) const;
private:
    pointer operator->() const;
};

typedef IteratorHelper<true> VertexIterator;
typedef IteratorHelper<false> EdgeIterator;

// implementation of Edge

template<bool DIRECTED, class T>
inline
Edge<DIRECTED, T>::Edge(
    const Value v0,
    const Value v1
) {
    if(DIRECTED) { // evaluated at compile time
        vertexIndices_[0] = v0;
        vertexIndices_[1] = v1;
    }
    else {
        if(v0 <= v1) {
            vertexIndices_[0] = v0;
            vertexIndices_[1] = v1;
        }
        else {
            vertexIndices_[0] = v1;
            vertexIndices_[1] = v0;
        }
    }
}

template<bool DIRECTED, class T>
inline typename Edge<DIRECTED, T>::Value
Edge<DIRECTED, T>::operator[]
(
    const std::size_t j
) const {
    assert(j < 2);

    return vertexIndices_[j];
}

template<bool DIRECTED, class T>
inline typename Edge<DIRECTED, T>::Value&
Edge<DIRECTED, T>::operator[](
    const std::size_t j
) {
    assert(j < 2);

    return vertexIndices_[j];
}

// implementation of IteratorHelper

template<bool T>
inline
IteratorHelper<T>::IteratorHelper()
:   Base()
{}

template<bool T>
inline
IteratorHelper<T>::IteratorHelper(
    const Base& it
)
:   Base(it)
{}

template<bool T>
inline
IteratorHelper<T>::IteratorHelper(
    const IteratorHelper<T>& it
)
:   Base(it)
{}

template<bool T>
inline IteratorHelper<T>
IteratorHelper<T>::operator=(
    const Base& it
) {
    Base::operator=(it);
    return *this;
}

template<bool T>
inline IteratorHelper<T>
IteratorHelper<T>::operator=(
    const IteratorHelper<T>& it
) {
    Base::operator=(it);
    return *this;
}

template<bool T>
inline typename IteratorHelper<T>::value_type
IteratorHelper<T>::operator*() const {
    if(T) { // evaluated at compile time
        return Base::operator*().vertex();
    }
    else {
        return Base::operator*().edge();
    }
}

template<bool T>
inline typename IteratorHelper<T>::value_type
IteratorHelper<T>::operator[](
    const std::size_t j
) const {
    if(T) { // evaluated at compile time
        return Base::operator[](j).vertex();
    }
    else {
        return Base::operator[](j).edge();
    }
}

template<bool T>
inline IteratorHelper<T>&
IteratorHelper<T>::operator+=(
    const difference_type d
) {
    Base::operator+=(d);
    return *this;
}

template<bool T>
inline IteratorHelper<T>&
IteratorHelper<T>::operator-=(
    const difference_type d
) {
    Base::operator-=(d);
    return *this;
}

template<bool T>
inline IteratorHelper<T>&
IteratorHelper<T>::operator++() { // prefix
    Base::operator++();
    return *this;
}

template<bool T>
inline IteratorHelper<T>&
IteratorHelper<T>::operator--() { // prefix
    Base::operator--();
    return *this;
}

template<bool T>
inline IteratorHelper<T>
IteratorHelper<T>::operator++(int) { // postfix
    return Base::operator++(int());
}

template<bool T>
inline IteratorHelper<T>
IteratorHelper<T>::operator--(int) { // postfix
    return Base::operator--(int());
}

template<bool T>
inline IteratorHelper<T>
IteratorHelper<T>::operator+(
    const difference_type d
) const {
    return Base::operator+(d);
}

template<bool T>
inline IteratorHelper<T>
IteratorHelper<T>::operator-(
    const difference_type d
) const {
    return Base::operator-(d);
}

#ifdef _MSC_VER
template<bool T>
inline typename IteratorHelper<T>::difference_type
IteratorHelper<T>::operator-(
    const IteratorHelper<T>& other
) const {
    return Base::operator-(other);
}
#endif

} // namespace detail
// \endcond
} // namespace graph
} // namespace andres

#endif // #ifndef ANDRES_GRAPH_DETAIL_GRAPH_HXX
