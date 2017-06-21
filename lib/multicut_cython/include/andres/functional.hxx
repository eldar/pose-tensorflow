#pragma once
#ifndef ANDRES_FUNCTIONAL_HXX
#define ANDRES_FUNCTIONAL_HXX

#include <cassert>
#include <cmath>
#include <stdexcept>

namespace andres {

template<typename T>
struct Identity {
    typedef T argument_type;
    typedef T result_type;
#ifndef _MSC_VER
	constexpr 
#endif
    result_type operator()(const argument_type& x) const
        { return x; }
};

template<typename T, typename U>
struct NegativeLogProbabilityRatio {
    typedef T argument_type;
    typedef U result_type;

    NegativeLogProbabilityRatio(const argument_type epsilon = static_cast<argument_type>(1) / static_cast<argument_type>(255))
        :   epsilon_(epsilon),
            oneMinusEpsilon_(static_cast<argument_type>(1) - epsilon)
        {
            if(epsilon <= 0 || epsilon * 2 >= 1) {
                throw std::out_of_range("epsilon out of range (0, 0.5).");
            }
        }
    result_type operator()(argument_type x) const
        {
            assert(0 <= x && x <= 1);
            if(x < epsilon_) {
                x = epsilon_;
            }
            else if(x > oneMinusEpsilon_) {
                x = oneMinusEpsilon_;
            }
            return std::log( (1-x)/x );
        }

private:
    const argument_type epsilon_;
    const argument_type oneMinusEpsilon_;
};

template<typename T, typename U>
struct NegativeLogProbabilityToInverseProbability {
    typedef T argument_type;
    typedef U result_type;

    result_type operator()(argument_type x) const {
        return 1-::exp(-x );
    }
};

template<typename T, typename U>
struct ProbabilityToNegativeLogInverseProbability {
    typedef T argument_type;
    typedef U result_type;

    result_type operator()(argument_type x) const {
        return -::log( 1-x );
    }
};

template<typename T, typename U>
struct ProbabilityToLogit {
    typedef T argument_type;
    typedef U result_type;

    result_type operator()(argument_type x) const {
      return ::log( (1-x)/x );
    }
};

} // namespace andres

#endif // #ifndef ANDRES_FUNCTIONAL_HXX
