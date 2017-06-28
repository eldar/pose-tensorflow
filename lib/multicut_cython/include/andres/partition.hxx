#pragma once
#ifndef ANDRES_PARTITION_HXX
#define ANDRES_PARTITION_HXX

#include <cstddef>
#include <vector>
#include <map>

/// The public API.
namespace andres {

/// Disjoint set data structure with path compression.
template<class T = std::size_t>
class Partition {
public:
    typedef T Index;

    Partition(const Index = 0);
    void assign(const Index = 0);

    Index find(const Index) const; // without path compression
    Index find(Index); // with path compression
    Index numberOfElements() const;
    Index numberOfSets() const;
    template<class Iterator>
        void elementLabeling(Iterator) const;
    template<class Iterator>
        void representatives(Iterator) const;
    void representativeLabeling(std::map<Index, Index>&) const;

    void merge(Index, Index);
    void insert(const Index);

private:
    std::vector<Index> parents_;
    std::vector<Index> ranks_;
    Index numberOfSets_;
};

/// Construct a partition (with a number of sets each containing one element).
///
/// \param size Number of distinct sets. 
///
template<class T>
inline 
Partition<T>::Partition(
    const Index size
)
:   parents_(static_cast<std::size_t>(size)),
    ranks_(static_cast<std::size_t>(size)),
    numberOfSets_(size)
{
    for(Index j = 0; j < size; ++j) {
        parents_[static_cast<std::size_t>(j)] = j;
    }
}

/// Reset the partition (to a number of sets each containing one element).
///
/// \param size Number of distinct sets.
///
template<class T>
inline void
Partition<T>::assign(
    const Index size
) {
    parents_.resize(static_cast<std::size_t>(size));
    ranks_.resize(static_cast<std::size_t>(size));
    numberOfSets_ = size;
    for(Index j = 0; j < size; ++j) {
        parents_[static_cast<std::size_t>(j)] = j;
    }
}

template<class T>
inline typename Partition<T>::Index 
Partition<T>::numberOfElements() const {
    return static_cast<Index>(parents_.size());
}

template<class T>
inline typename Partition<T>::Index 
Partition<T>::numberOfSets() const {
    return numberOfSets_; 
}

/// Find the representative element of the set that contains the given element (without path compression).
/// 
/// \param element Element. 
///
template<class T>
inline typename Partition<T>::Index
Partition<T>::find(
    const Index element
) const {
    // find the root
    Index root = element;
    while(parents_[static_cast<std::size_t>(root)] != root) {
        root = parents_[static_cast<std::size_t>(root)];
    }
    return root;
}

/// Find the representative element of the set that contains the given element (with path compression).
/// 
/// This mutable function compresses the search path.
///
/// \param element Element. 
///
template<class T>
inline typename Partition<T>::Index
Partition<T>::find(
    Index element // copy to work with
) {
    // find the root
    Index root = element;
    while(parents_[static_cast<std::size_t>(root)] != root) {
        root = parents_[static_cast<std::size_t>(root)];
    }
    // path compression
    while(element != root) {
        const Index tmp = parents_[static_cast<std::size_t>(element)];
        parents_[static_cast<std::size_t>(element)] = root;
        element = tmp;
    }
    return root;
}

/// Merge two sets.
/// 
/// \param element1 Element in the first set. 
/// \param element2 Element in the second set. 
///
template<class T>
inline void 
Partition<T>::merge(
    Index element1,
    Index element2
) {
    // merge by rank
    element1 = find(element1);
    element2 = find(element2);
    if(ranks_[static_cast<std::size_t>(element1)] < ranks_[static_cast<std::size_t>(element2)]) {
        parents_[static_cast<std::size_t>(element1)] = element2;
        --numberOfSets_;
    }
    else if(ranks_[static_cast<std::size_t>(element1)] > ranks_[static_cast<std::size_t>(element2)]) {
        parents_[static_cast<std::size_t>(element2)] = element1;
        --numberOfSets_;
    }
    else if(element1 != element2) {
        parents_[static_cast<std::size_t>(element2)] = element1;
        ++ranks_[static_cast<std::size_t>(element1)];
        --numberOfSets_;
    }
}

/// Insert a number of new sets, each containing one element.
/// 
/// \param number Number of sets to insert. 
///
template<class T>
inline void 
Partition<T>::insert(
    const Index number
) {
    const Index numberOfElements = static_cast<Index>(parents_.size());
    ranks_.insert(ranks_.end(), static_cast<std::size_t>(number), 0);
    parents_.insert(parents_.end(), static_cast<std::size_t>(number), 0);
    for(Index j = numberOfElements; j < numberOfElements + number; ++j) {
        parents_[static_cast<std::size_t>(j)] = j;
    }
    numberOfSets_ += number;
}

/// Output all elements which are set representatives.
/// 
/// \param it (Output) Iterator.
///
template<class T>
template<class Iterator>
inline void 
Partition<T>::representatives(
    Iterator it
) const {
    for(Index j = 0; j < numberOfElements(); ++j) {
        if(parents_[static_cast<std::size_t>(j)] == j) {
            *it = j;
            ++it;
        }
    }
}

/// Output a contiguous labeling of the representative elements.
/// 
/// \param out (Output) A map that assigns each representative element to an integer label.
///
template<class T>
inline void 
Partition<T>::representativeLabeling(
    std::map<Index, Index>& out
) const {
    out.clear();	
    std::vector<Index> r(static_cast<std::size_t>(numberOfSets()));
    representatives(r.begin());
    for(Index j = 0; j < numberOfSets(); ++j) {
        out[r[static_cast<std::size_t>(j)]] = j;
    }
}

/// Output a contiguous labeling of all elements.
/// 
/// \param out (Output) Iterator into a container in which the j-th entry becomes the label of the j-th element.
///
template<class T>
template<class Iterator>
inline void 
Partition<T>::elementLabeling(
    Iterator out
) const {
    std::map<Index, Index> rl;
    representativeLabeling(rl);
    for(Index j = 0; j < numberOfElements(); ++j) {
        *out = rl[find(j)];
        ++out;
    }
}

} // namespace andres

#endif // #ifndef ANDRES_PARTITION_HXX
