/*
 * Taken from http://www.andrew.cmu.edu/user/vanhoeve/mdd/ with the notice:
 * "The software can be freely used but comes with no warranty"
 *
 * bdd.hpp
 *      Author: Bergman, Cire
 */

#ifndef BDD22_HPP_
#define BDD22_HPP_

#include <map>
#include <vector>
#include "intset.hpp"
#include "util.hpp"

using namespace std;



/**
 * Node2
 */
struct Node2 {

//	int				id;

	IntSet			state;
	int				longest_path;
//	bool			exact;

	int				relax_ub;


	/**
	 * Node2 constructor if one wishes only to create a relaxation
	 */
	Node2(IntSet &_state, int _longest_path)	: state(_state), longest_path(_longest_path)
	{
	};
};	


struct BDD2 {
};


typedef map<IntSet*, Node2*, IntSetLexLessThan> Node2Map2;


/**
 * Node2 comparator by longest path
 */
struct CompareNode2sLongestPath {
	bool operator()(const Node2* nodeA, const Node2* nodeB) const {
		return nodeA->longest_path > nodeB->longest_path;
	}
};

/**
 * Node2 comparator by state size
 */
struct CompareNode2sStateSizeAscending {
	bool operator()(Node2* nodeA, Node2* nodeB) const {
		if( nodeA->state.get_size() != nodeB->state.get_size() )
			return nodeA->state.get_size() < nodeB->state.get_size();
		return nodeA->longest_path > nodeB->longest_path;
	}
};

/**
 * Node2 comparator by state size
 */
struct CompareNode2sStateSizeDescending {
	bool operator()(Node2* nodeA, Node2* nodeB) const {
		if( nodeA->state.get_size() != nodeB->state.get_size() )
			return nodeA->state.get_size() > nodeB->state.get_size();
		return nodeA->longest_path > nodeB->longest_path;
	}
};



#endif /* BDD2_HPP_ */
