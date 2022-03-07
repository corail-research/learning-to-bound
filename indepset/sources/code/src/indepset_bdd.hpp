/*
 * ===============================================================
 * IndepSetBDD - Definition class
 * ===============================================================
 */

#ifndef INDEPSETBDD_HPP_
#define INDEPSETBDD_HPP_

#include <map>
#include <queue>
#include <deque>
#include <string>
#include <unistd.h>

#include "indepset_orderings.hpp"
#include "indepset_instance.hpp"
#include "util.hpp"

#include "config.h"
#include "learning_lib.h"
#include "graph.h"
#include "nn_api.h"
#include "misp_qnet.h"
#include "nstep_replay_mem.h"
#include "simulator.h"
#include "learning_env.h"

#define EPS 1e-8

using namespace std;


//
// Parameters
// TODO: include this in a parameter struct
//
#define ROOT_WIDTH	100				// width for initial relaxation

//extern GSet GSetTest2;

// -------------------------------------------------------------------
// Branching Node
// -------------------------------------------------------------------

struct BranchNode {
	vector<int>		state;				// state to be explored
	int   			longest_path;		// node longest path
	int   			relax_ub;           // upper bound obtained for the BDD where node was created
	int			oracle_added;		//added by oracle if == 1

	//
	// Constructor from a BDD node
	//
	BranchNode(Node* _node, double _longest_path, int _o_a)
	:	longest_path(_longest_path), relax_ub(0), oracle_added(_o_a)
	{
		const State& state_v = _node->state;
		state.reserve( state_v.count() );
		for (int v = state_v.find_first(); v != state_end; v = state_v.find_next(v)) {
			state.push_back(v);
		}
	}

	//
	// Constructor from a state vector
	//
	BranchNode(const vector<int>& _state, double _longest_path)
	:	state(_state), longest_path(_longest_path), relax_ub(0)
	{
	}

	//
	// Empty constructor
	//
	BranchNode() {
	}
};


//
// Branch node comparator by relaxed UB
// b1 comes before b2 if relaxUB is lower 
//
struct BranchNodeComparatorUB {
	BranchNodeComparatorUB() { }
	bool operator()(const BranchNode* b1, const BranchNode* b2) {
		if (b1->relax_ub == b2->relax_ub) {
			if (b1->state.size() == b2->state.size()) {
				return b1->longest_path < b2->longest_path;
			}
			return b1->state.size() < b2->state.size();
		}
		return b1->relax_ub < b2->relax_ub;
	}
};

//
// Branch node comparator for DFS
//
struct BranchNodeComparatorDFS {
	BranchNodeComparatorDFS() { }
	bool operator()(const BranchNode* b1, const BranchNode* b2) {
		if (b1->state.size() == b2->state.size()) {
			if (b1->relax_ub == b2->relax_ub) {
				return b1->longest_path < b2->longest_path;
			}
			return b1->relax_ub < b2->relax_ub;
		}
		return b1->state.size() > b2->state.size();
	}
};


//
// Queue of branch nodes ordered by relax ub
//
typedef priority_queue<BranchNode*, vector<BranchNode*>, BranchNodeComparatorUB> BranchNodeQueue;


//
// Queue of branch nodes ordered by state size for DFS search
//
typedef priority_queue<BranchNode*, vector<BranchNode*>, BranchNodeComparatorDFS> BranchNodeDFS;


// -------------------------------------------------------------------
// Independent Set BDD: Relaxation / Restriction generator
// -------------------------------------------------------------------

// Independent set solver class
class IndepSetBDD {

public:

	//
	// Constructor given a DIMACS instance filename
	//
	// IndepSetBDD(char* filename, int _maxWidth)
	// DDX10:
	IndepSetBDD(const int _rootWidth,
				const int _ddWidth,
				const string & _instanceName,
				std::string filename,
				const int argc,
				char** args,
				int oracle,
				int rand,
				int nqval);

	// Generate relaxation. Returns true iff relaxation is exact.
	int generateRelaxationOracle(const int initial_weight, std::string model, const int argc, char** args, int nqval);

	int generateRelaxationHeuristic(int initial_lp,int rand);

	// Get ordering pre-computed depending on the depth and variables covered
	void GetOrder(int depth, int *sol, vector<int> act_ver, std::string model, const int argc, char** args);

	// Generate relaxation given a Branch node
	int generateRelaxation(const BranchNode* branch_node, std::string model, const int argc, char** args, int oracle, int rand, int nqval);

	// Generate restriction
	int generateRestriction(const int initial_weight);

	// Generate restriction given a Branch node
	int generateRestriction(const BranchNode* branch_node);

	// Update local lower bound, indicating whether it was found locally or remotely
	void updateLocalLB(int _lb, bool isLocal = false);

	// Verify if relaxed BDD is exact
	bool isBDDExact() 				{ return isExact; }

	// Get best lower bound found
	int getBestLB() 		  		{ return best_lb; }

	// Get last computed upper bound
	int getUB() 		  			{ return upper_bound; }

	// Get size
	int getSizeVert() 		  			{ return n_var ; }
	int getSizeBranch() 		  			{ return branch_nodes.size() ; }

	// Add branching nodes to queue
  void addBranchNodesQueue(BranchNodeQueue& queue, unsigned long int &size_pool);

	// Add branching nodes to DFS queue
	void addBranchNodesDFS(BranchNodeDFS& queue);

	// Add branching nodes to DFS queue
  void addBranchNodesDFS(vector<BranchNode*>& queue);


	// Public data members
	IndepSetInst* 			Inst;				// independent set instance
	vector<int> 			vertex_in_layer;  	// which vertex a layer corresponds to
	Node*					root_node;		  	// relaxation root node

	int n_var;

	int ask_oracle_count = 0;

	// for pre-computed orderings, d and orderings
	deque<int> d_orderings; //depths of pre-computed orderings
	deque<vector<int>> orderings; //corresponding orderings
	int maxLenOrd = 500;

private:

	// ----------------------------------------------
	// Attributes
	// ----------------------------------------------

	BDDNodePool 		node_pool;		  		// list of nodes to explore

	vector<int>			active_vertices;  		// active vertices in the relaxation
	vector<Node*> 		nodes_layer;	  		// nodes in a layer

	int* 				in_state_counter;  		// min-in-state vertex ordering parameters
	IS_Ordering* 		var_ordering; 			// BDD Vertex ordering

	int 				maxWidth;  				// max allowed width
	double 				best_lb;				// best lower bound
	bool				isExact;				// if BDD is exact

	bool				isLBUpdated;			// if lower bound has been updated

	vector<BranchNode*> branch_nodes; 			// BDD nodes that require branching

	vector<int>			root_ordering;			// relative ordering of each vertex

	int					upper_bound;			// last computed upper bound

	int					initialPosPool;

	
	

	// Auxiliary parameters
	CompareNodesLongestPath bn_comparator;


	// ----------------------------------------------
	// Functions
	// ----------------------------------------------

	// Choose vertex
	int choose_next_vertex_min_size_next_layer();

	// Merge two nodes
	void merge(Node* nodeA, Node* nodeB, int oa);

	// Merge layer
	void mergeLayer(int layer, vector<Node*> &nodes_layer, int oa);

	// Restrict layer
	void restrictLayer(int layer, vector<Node*> &nodes_layer);

	// Add new local node to be explored later
	void addBranchNode(Node* _node, int oa);
};


// -------------- Utilities ----------------

struct CompareOrderDescending {
	vector<int> &order;
	CompareOrderDescending(vector<int>& _order) : order(_order) { }
	bool operator()(const int a, const int b) {
		return order[a] > order[b];
	}
};


// ----------------------------------------------------
// Inline implementations
// ----------------------------------------------------

//
// Merge node A into node B. Node A is always deleted.
//
inline void IndepSetBDD::merge(Node* nodeA, Node* nodeB, int oa) {
	assert( (nodeA->state) == (nodeB->state) );

	if ( (nodeA->exact == nodeB->exact) ) {

		// Both nodes are either relaxed or exact
		nodeB->longest_path = MAX(nodeA->longest_path, nodeB->longest_path);

	} else if (nodeA->exact) {

		// Node A is exact, B is relaxed
		branch_nodes.push_back( new BranchNode(nodeA, nodeA->longest_path, oa) );
		nodeB->longest_path = MAX(nodeA->longest_path, nodeB->longest_path);

	} else if (nodeB->exact) {

		// Node B is exact, A is relaxed
		branch_nodes.push_back( new BranchNode(nodeB, nodeB->longest_path, oa));

		nodeB->exact = false;
		nodeB->longest_path = MAX(nodeA->longest_path, nodeB->longest_path);
	}

	delete nodeA;
}

//
// Update lower bound, indicating whether it was found locally or remotely. If locally,
// then mark that bound was updated to send it to other nodes
//
inline void IndepSetBDD::updateLocalLB(int _lb, bool isLocal) {
	if (_lb > best_lb) {
		best_lb = _lb;
		if (isLocal) {
			isLBUpdated = true;
		}
	}
}


//
// Add new local node to be explored
//
inline void IndepSetBDD::addBranchNode(Node* _node, int oa) {
	branch_nodes.push_back(	new BranchNode(_node, _node->longest_path, oa) );
}

//
// Add branching nodes to queue
// TODO: apply bounding heuristic
//
inline void IndepSetBDD::addBranchNodesQueue(BranchNodeQueue& queue, unsigned long int &size_pool) {
	for (vector<BranchNode*>::iterator st = branch_nodes.begin(); st != branch_nodes.end(); ++st) {
		if ( (*st)->relax_ub <= best_lb ) {
			delete (*st);
		} else {
		//cout << "+1" << endl;
		  size_pool += (*st)->state.size();
		  queue.push(*st);
		}
	}
	branch_nodes.clear();
}


//
// Add branching nodes to DFS queue
// TODO: apply bounding heuristic
//
inline void IndepSetBDD::addBranchNodesDFS(BranchNodeDFS& queue) {
	for (vector<BranchNode*>::iterator st = branch_nodes.begin(); st != branch_nodes.end(); ++st) {
		if ( (*st)->relax_ub <= best_lb ) {
			delete (*st);
		} else {
			queue.push(*st);
		}
	}
	branch_nodes.clear();
}

//
// Add branching nodes to DFS queue
// TODO: apply bounding heuristic
//
inline void IndepSetBDD::addBranchNodesDFS(vector<BranchNode*>& queue) {
	for (vector<BranchNode*>::iterator st = branch_nodes.begin(); st != branch_nodes.end(); ++st) {
		if ( (*st)->relax_ub <= best_lb ) {
			delete (*st);
		} else {
			queue.push_back(*st);
		}
	}
	branch_nodes.clear();
}


//
// Generate relaxation given a branch node
//
inline int IndepSetBDD::generateRelaxation(const BranchNode* branch_node, std::string model, const int argc, char** args, int oracle, int rand, int nqval) {

	active_vertices = branch_node->state;
	if(oracle){
	return generateRelaxationOracle(branch_node->longest_path, model, argc, args, nqval);}
	else{return generateRelaxationHeuristic(branch_node->longest_path, rand);}
}

//
// Generate restriction given a branch node
//
inline int IndepSetBDD::generateRestriction(const BranchNode* branch_node) {
	active_vertices = branch_node->state;
	return generateRestriction(branch_node->longest_path);
}


#endif /* INDEPSETBDD_HPP_ */
