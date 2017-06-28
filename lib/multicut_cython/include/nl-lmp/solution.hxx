#pragma once
#ifndef NL_LMP_SOLUTION_HXX
#define NL_LMP_SOLUTION_HXX

#include <vector>


namespace nl_lmp
{

class Solution
{
public:
    struct element
    {
        size_t classIndex { 0 };
        size_t clusterIndex { 0 };
    };

    explicit Solution(size_t n) :
        solution_(n)
    {}

    Solution(Solution const& other)
    {
        solution_ = other.solution_;
    }

    element const& operator[](size_t index) const
    {
        return solution_[index];
    }

    std::vector<element>::const_iterator begin() const
    {
        return solution_.cbegin();
    }

    std::vector<element>::const_iterator end() const
    {
        return solution_.cend();
    }

    size_t size() const
    {
        return solution_.size();
    }

    element& operator[](size_t index)
    {
        return solution_[index];
    }

    std::vector<element>::iterator begin()
    {
        return solution_.begin();
    }

    std::vector<element>::iterator end()
    {
        return solution_.end();
    }

private:
    std::vector<element> solution_;
};

}

#endif
