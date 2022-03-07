/*
 * ------------------------------------
 * General utility functions
 * ------------------------------------
 */

#ifndef UTIL_HPP_
#define UTIL_HPP_

#include <limits>


/**
 * -------------------------------------------------------------
 * Macros
 * -------------------------------------------------------------
 */

#define MIN(_a_,_b_)              ((_a_ < _b_) ? _a_ : _b_)
#define MAX(_a_,_b_)              ((_a_ > _b_) ? _a_ : _b_)
#define FLOOR(_a_)                ( (int)_a_ )
#define CEIL(_a_)                 ( ((int)_a_ == _a_) ? _a_ : (((int)_a_)+1) )
//#define ABS(_a_)                    ( _a_ < 0 ? (-1)*_a_ : _a_ )


/**
 * -------------------------------------------------------------
 * Constants
 * -------------------------------------------------------------
 */

const int   INF_WIDTH   = -1;        							 // infinite width
const int   INF   		= std::numeric_limits<int>::max();		 // infinite *


/**
 * -------------------------------------------------------------
 * Comparator with respect to auxiliary vector
 * -------------------------------------------------------------
 */

struct ComparatorAuxIntVectorDescending {
	const std::vector<int> &v;
	ComparatorAuxIntVectorDescending(const std::vector<int> &_v) : v(_v) { }
	bool operator()(int i, int j) {
		return (v[i] > v[j]);
	}
};


struct ComparatorAuxIntVectorAscending {
	const std::vector<int> &v;
	ComparatorAuxIntVectorAscending(const std::vector<int> &_v) : v(_v) { }
	bool operator()(int i, int j) {
		return (v[i] > v[j]);
	}
};


#endif /* UTIL_HPP_ */
