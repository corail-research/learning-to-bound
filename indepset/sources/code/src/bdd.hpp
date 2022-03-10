/*
 * ---------------------------------------------------------
 * BDD Structure definitions for independet set
 * ---------------------------------------------------------
 */

#ifndef BDD_HPP_
#define BDD_HPP_


#include <map>
#include <vector>
#include <boost/dynamic_bitset.hpp>
#include "util.hpp"

using namespace std;


typedef boost::dynamic_bitset<> State;


//
// Lexicographic state comparator
//
struct StateLessThan {
	bool operator()(const State* stateA, const State* stateB) const {
		return (*stateA) < (*stateB);
	}
};


//
// BDD Node
//
struct Node {
	State			state;
	double			longest_path;
	bool			exact;

	Node*			one_arc;
	Node*			zero_arc;
	int 			layer;


	Node(State &_state, double _longest_path)
	: state(_state), longest_path(_longest_path), exact(false), one_arc(NULL), zero_arc(NULL), layer(-1)
	{ }

	Node(State &_state, double _longest_path, bool _exact)
	: state(_state), longest_path(_longest_path), exact(_exact), one_arc(NULL), zero_arc(NULL), layer(-1)
	{ }
};

//
// BDD data structure
//
struct BDD {

};

//
// BDD Node pool type
//
typedef map<State*, Node*, StateLessThan> BDDNodePool;

//
// Marker of the end of a state (for iteration purposes)
//
const int state_end = static_cast<int>(boost::dynamic_bitset<>::npos);


/*
 * ---------------------------------------------------------
 * BDD Node comparators
 * ---------------------------------------------------------
 */

//
// Node comparator by longest path
//
struct CompareNodesLongestPath {
	bool operator()(const Node* nodeA, const Node* nodeB) const {
		return nodeA->longest_path > nodeB->longest_path;
	}
};


#endif /* BDD_HPP_ */
