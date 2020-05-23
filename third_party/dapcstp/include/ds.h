/**
 * \file   ds.h
 * \brief  data structure definitions
 *
 * \author Martin Luipersbeck 
 * \date   2015-05-03
 */

#ifndef DS_H_
#define DS_H_

#include "def.h"

#include <boost/heap/d_ary_heap.hpp>
#include <boost/heap/fibonacci_heap.hpp>

// stores information about transformation from MWCS to PCSTP
struct Transformation {
	weight_t P, minP, sumC;
};

struct InstSizeData {
	int n, m, t, tr, f1;
};

// comparison operators (priority queues return top (last) element)
template<typename T, typename U>
struct largestFirst {
	bool operator() (const pair<T,U> a, const pair<T,U> b) const {
		return a.first < b.first;
	}
};

template<typename T, typename U>
struct smallestFirst {
	bool operator() (const pair<T,U> a, const pair<T,U> b) const {
		return a.first > b.first;
	}
};
template <typename T, typename U> using PQMin = boost::heap::d_ary_heap<pair<T,U>, boost::heap::arity<2>, boost::heap::compare<smallestFirst<T,U>>, boost::heap::mutable_<true>>;

template <typename T, typename U> using PQMax = boost::heap::d_ary_heap<pair<T,U>, boost::heap::arity<2>, boost::heap::compare<largestFirst<T,U>>, boost::heap::mutable_<true>>;


template <typename T, typename U> using PQMinFib = boost::heap::fibonacci_heap<pair<T,U>, boost::heap::mutable_<true>, boost::heap::compare<smallestFirst<T,U>>>;

#endif // DS_H_
