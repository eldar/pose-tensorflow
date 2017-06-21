#pragma once
#ifndef ANDRES_RANDOM_ACCESS_SET_HXX
#define ANDRES_RANDOM_ACCESS_SET_HXX

#include <cstddef>
#include <vector>
#include <algorithm> // std::lower_bound, std::upper_bound
#include <functional> // std::greater
#include <utility> // std::make_pair

namespace andres {

// STL-compliant container compatible with std::set but with linear time insert and constant time random access.
//
// this file is based on randomaccessset.hxx written by Thorsten Beier for OpenGM: 
// http://hci.iwr.uni-heidelberg.de/opengm2/
//
template<
    class Key, 
    class Comparison = std::less<Key>, 
    class Allocator = std::allocator<Key> 
>
class RandomAccessSet {
private:
    typedef std::vector<Key, Allocator> Vector;

public:
    typedef Key key_type;
    typedef Key value_type;
    typedef Comparison key_compare;
    typedef Comparison value_compare;
    typedef Allocator allocator_type;

    typedef typename Allocator::const_reference const_reference;
    typedef typename Vector::iterator iterator;
    typedef typename Vector::const_iterator const_iterator;
    typedef typename Vector::reverse_iterator reverse_iterator;
    typedef typename Vector::const_reverse_iterator const_reverse_iterator;
    typedef typename Vector::size_type size_type;
    typedef typename Vector::difference_type	difference_type;
    typedef typename Vector::const_pointer const_pointer;

    RandomAccessSet(const std::size_t, const Comparison& = Comparison(),
        const Allocator& = Allocator());
    RandomAccessSet(const Comparison& = Comparison(), 
        const Allocator& = Allocator());
    template <class Iterator>
    RandomAccessSet(Iterator, Iterator, const Comparison& = Comparison(), 
        const Allocator& = Allocator());

    const value_type& operator[](const size_type) const;

    const_iterator begin() const;
    const_iterator end() const;
    iterator begin();
    iterator end();

    const_reverse_iterator rbegin() const;
    const_reverse_iterator rend() const;
    reverse_iterator rbegin();
    reverse_iterator rend();

    const_iterator find(const key_type&) const;
    iterator find(const key_type&);

    bool empty() const;
    size_type size() const;
    size_type max_size() const;

    std::pair<const_iterator, bool> insert(const value_type&);
    template <class Iterator>
    void insert(Iterator, Iterator);
    const_iterator insert(iterator, const value_type&); 

    void erase(iterator position);
    size_type erase(const key_type&);
    void erase(iterator, iterator);   
    void clear();

    size_type count(const key_type&) const;
    key_compare key_comp() const;
    value_compare value_comp() const;
    const_iterator lower_bound(const key_type&) const;
    const_iterator upper_bound(const key_type&) const;
    iterator lower_bound(const key_type&);
    iterator upper_bound(const key_type&);
    std::pair<const_iterator, const_iterator> equal_range(const key_type&) const;
    std::pair<iterator, iterator> equal_range(const key_type&);

    allocator_type get_allocator() const;

    // TODO: implement C++11 member functions 'emplace' and 'emplace_hint'

private:
    std::vector<Key> vector_;
    Comparison compare_;
};

template<class Key, class Comparison, class Allocator>
inline
    RandomAccessSet<Key, Comparison, Allocator>::RandomAccessSet(
    const Comparison& comparison, 
    const Allocator& allocator
)
:   vector_(allocator), 
    compare_(comparison)
{}

template<class Key, class Comparison, class Allocator>
inline
    RandomAccessSet<Key, Comparison, Allocator>::RandomAccessSet(
    const std::size_t reserveSize,
    const Comparison& comparison, 
    const Allocator& allocator
)
:   vector_(allocator), 
    compare_(comparison) 
{
    vector_.reserve(reserveSize);
}

template<class Key, class Comparison, class Allocator>
template <class Iterator>
inline
    RandomAccessSet<Key, Comparison, Allocator>::RandomAccessSet(
    Iterator beginInput, 
    Iterator endInput, 
    const Comparison& comparison, 
    const Allocator& allocator
)
:   vector_(allocator), 
    compare_(comparison)
{
    while(beginInput != endInput) {
        insert(*beginInput);
        ++beginInput;
    }
}

template<class Key, class Comparison, class Allocator>
inline const typename RandomAccessSet<Key, Comparison, Allocator>::value_type&
    RandomAccessSet<Key, Comparison, Allocator>::operator[](
    const typename RandomAccessSet<Key, Comparison, Allocator>::size_type  index
) const {
    return vector_[index];
}

template<class Key, class Comparison, class Allocator>
inline typename RandomAccessSet<Key, Comparison, Allocator>::const_iterator
RandomAccessSet<Key, Comparison, Allocator>::begin() const {
    return vector_.begin();
}

template<class Key, class Comparison, class Allocator>
inline typename RandomAccessSet<Key, Comparison, Allocator>::const_iterator
RandomAccessSet<Key, Comparison, Allocator>::end() const {
    return vector_.end();
}

template<class Key, class Comparison, class Allocator>
inline typename RandomAccessSet<Key, Comparison, Allocator>::const_reverse_iterator
RandomAccessSet<Key, Comparison, Allocator>::rbegin() const {
    return vector_.rbegin();
}

template<class Key, class Comparison, class Allocator>
inline typename RandomAccessSet<Key, Comparison, Allocator>::const_reverse_iterator
RandomAccessSet<Key, Comparison, Allocator>::rend() const {
    return vector_.rend();
}

template<class Key, class Comparison, class Allocator>
inline typename RandomAccessSet<Key, Comparison, Allocator>::iterator
RandomAccessSet<Key, Comparison, Allocator>::begin() {
    return vector_.begin();
}

template<class Key, class Comparison, class Allocator>
inline typename RandomAccessSet<Key, Comparison, Allocator>::iterator
RandomAccessSet<Key, Comparison, Allocator>::end() {
    return vector_.end();
}

template<class Key, class Comparison, class Allocator>
inline typename RandomAccessSet<Key, Comparison, Allocator>::reverse_iterator
RandomAccessSet<Key, Comparison, Allocator>::rbegin() {
    return vector_.rbegin();
}

template<class Key, class Comparison, class Allocator>
inline typename RandomAccessSet<Key, Comparison, Allocator>::reverse_iterator
RandomAccessSet<Key, Comparison, Allocator>::rend() {
    return vector_.rend();
}

template<class Key, class Comparison, class Allocator>
inline bool
RandomAccessSet<Key, Comparison, Allocator>::empty() const {
    return vector_.empty();
}

template<class Key, class Comparison, class Allocator>
inline typename RandomAccessSet<Key, Comparison, Allocator>::size_type
RandomAccessSet<Key, Comparison, Allocator>::size() const {
    return vector_.size();
}

template<class Key, class Comparison, class Allocator>
inline typename RandomAccessSet<Key, Comparison, Allocator>::size_type
RandomAccessSet<Key, Comparison, Allocator>::max_size() const {
    return vector_.max_size();
}

template<class Key, class Comparison, class Allocator>
inline std::pair<typename RandomAccessSet<Key, Comparison, Allocator>::const_iterator, bool>
RandomAccessSet<Key, Comparison, Allocator>::insert(
    const typename RandomAccessSet<Key, Comparison, Allocator>::value_type& value
) {
    bool found(true);
    iterator i(lower_bound(static_cast<Key>(value)));
    if(i == end() || compare_(static_cast<Key>(value), *i)) {
        i = vector_.insert(i, static_cast<Key>(value));
        found = false;
    }
    return std::make_pair(i, !found);
}

template<class Key, class Comparison, class Allocator>
template <class Iterator>
inline void
    RandomAccessSet<Key, Comparison, Allocator>::insert(
    Iterator first, 
    Iterator last
) {
    for(; first != last; ++first) {
        insert(*first);
    }
}

// TODO: optimize according to C++11 specification:
// The function optimizes its insertion time if position points to the 
// element that will follow the inserted element (or to the end, if it 
// would be the last).
//
template<class Key, class Comparison, class Allocator>
inline typename RandomAccessSet<Key, Comparison, Allocator>::const_iterator
    RandomAccessSet<Key, Comparison, Allocator>::insert(
    typename RandomAccessSet<Key, Comparison, Allocator>::iterator position, 
    const typename RandomAccessSet<Key, Comparison, Allocator>::value_type& value
) {
    std::pair<const_iterator, bool> ret;
    ret = insert(value);
    return ret.first;
}

template<class Key, class Comparison, class Allocator>
inline void
RandomAccessSet<Key, Comparison, Allocator>::erase(
    typename RandomAccessSet<Key, Comparison, Allocator>::iterator position
) {
    vector_.erase(position);
}

template<class Key, class Comparison, class Allocator>
inline typename RandomAccessSet<Key, Comparison, Allocator>::size_type
RandomAccessSet<Key, Comparison, Allocator>::erase(
    const key_type& x
) {
    iterator i = find(x);
    if(i != vector_.end()) {
        erase(i);
        return 1;
    }
    return 0;
}

template<class Key, class Comparison, class Allocator>
inline void
RandomAccessSet<Key, Comparison, Allocator>::erase(
    const typename RandomAccessSet<Key, Comparison, Allocator>::iterator first, 
    const typename RandomAccessSet<Key, Comparison, Allocator>::iterator last
) {
    vector_.erase(first, last);
}

template<class Key, class Comparison, class Allocator>
inline void
RandomAccessSet<Key, Comparison, Allocator>::clear() {
    vector_.clear();
}

template<class Key, class Comparison, class Allocator>
inline typename RandomAccessSet<Key, Comparison, Allocator>::key_compare
RandomAccessSet<Key, Comparison, Allocator>::key_comp() const {
    return compare_;
}

template<class Key, class Comparison, class Allocator>
inline typename RandomAccessSet<Key, Comparison, Allocator>::value_compare
RandomAccessSet<Key, Comparison, Allocator>::value_comp() const {
    return compare_;
}

template<class Key, class Comparison, class Allocator>
inline typename RandomAccessSet<Key, Comparison, Allocator>::const_iterator
RandomAccessSet<Key, Comparison, Allocator>::find(
    const key_type& value
) const {
    const_iterator i(lower_bound(value));
    if(i != end() && compare_(value, *i)) {
        i = end();
    }
    return i;
}

template<class Key, class Comparison, class Allocator>
inline typename RandomAccessSet<Key, Comparison, Allocator>::iterator
RandomAccessSet<Key, Comparison, Allocator>::find(
    const key_type& value
) {
    iterator i(lower_bound(value));
    if(i != end() && compare_(value, *i)) {
        i = end();
    }
    return i;
}

template<class Key, class Comparison, class Allocator>
inline typename RandomAccessSet<Key, Comparison, Allocator>::size_type
RandomAccessSet<Key, Comparison, Allocator>::count(
    const key_type&  value
) const {
    if(find(value) == end()) {
        return 0;
    }
    else {
        return 1;
    }
}

template<class Key, class Comparison, class Allocator>
inline typename RandomAccessSet<Key, Comparison, Allocator>::const_iterator
RandomAccessSet<Key, Comparison, Allocator>::lower_bound(
    const key_type& value
) const {
    return std::lower_bound(vector_.begin(), vector_.end(), value, compare_);
}

template<class Key, class Comparison, class Allocator>
inline typename RandomAccessSet<Key, Comparison, Allocator>::iterator
RandomAccessSet<Key, Comparison, Allocator>::lower_bound(
    const key_type& value
) {
    return std::lower_bound(vector_.begin(), vector_.end(), value, compare_);
}

template<class Key, class Comparison, class Allocator>
inline typename RandomAccessSet<Key, Comparison, Allocator>::const_iterator
RandomAccessSet<Key, Comparison, Allocator>::upper_bound(
    const key_type& value
) const {
    return std::upper_bound(vector_.begin(), vector_.end(), value, compare_);
}

template<class Key, class Comparison, class Allocator>
inline typename RandomAccessSet<Key, Comparison, Allocator>::iterator
RandomAccessSet<Key, Comparison, Allocator>::upper_bound(
    const key_type& value
) {
    return std::upper_bound(vector_.begin(), vector_.end(), value, compare_);
}

template<class Key, class Comparison, class Allocator>
inline std::pair<typename RandomAccessSet<Key, Comparison, Allocator>::const_iterator, typename RandomAccessSet<Key, Comparison, Allocator>::const_iterator>
RandomAccessSet<Key, Comparison, Allocator>::equal_range(
    const key_type& value
) const {
    return std::equal_range(vector_.begin(), vector_.end(), value, compare_);
}

template<class Key, class Comparison, class Allocator>
inline std::pair<typename RandomAccessSet<Key, Comparison, Allocator>::iterator, typename RandomAccessSet<Key, Comparison, Allocator>::iterator>
RandomAccessSet<Key, Comparison, Allocator>::equal_range(
    const key_type& value
) {
    return std::equal_range(vector_.begin(), vector_.end(), value, compare_);
}

template<class Key, class Comparison, class Allocator>
inline typename RandomAccessSet<Key, Comparison, Allocator>::allocator_type
RandomAccessSet<Key, Comparison, Allocator>::get_allocator() const {
    return vector_.get_allocator();
}

} // namespace andres

#endif // #ifndef ANDRES_RANDOM_ACCESS_SET_HXX
