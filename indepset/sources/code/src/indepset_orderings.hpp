/*
 * --------------------------------------------------------
 * BDD Vertex Ordering class
 * --------------------------------------------------------
 */

#ifndef ORDERING_HPP_
#define ORDERING_HPP_

#include <cassert>
#include <cstdio>
#include <vector>

#include "indepset_instance.hpp"
#include "bdd.hpp"

using namespace std;


struct IntComparator {
	vector<int> &v;
	IntComparator(vector<int> &_v) : v(_v) { }
	bool operator()(int i, int j) {
		return (v[i] > v[j]);
	}
};


using namespace std;

enum OrderType {
  MinState, MaximalPath, CutVertexGen, CutVertex, Fixed,
  MinDegree, RootOrder, LexOrder, RandMinState, Random
};

// Class representing a general ordering
struct IS_Ordering {

	IndepSetInst   *inst;
	char           name[256];
	OrderType	   order_type;

	IS_Ordering(IndepSetInst* _inst, OrderType _order_type) : inst(_inst), order_type(_order_type) { }

	virtual ~IS_Ordering() { }

	// returns vertex corresponding to particular layer
	virtual int vertex_in_layer(int layer) = 0;
};



// Ordering according to root node
struct MinInState : IS_Ordering {

	vector<bool> avail_v;   // available vertices

	MinInState(IndepSetInst *_inst) : IS_Ordering(_inst, MinState) {
		avail_v.resize(inst->graph->n_vertices, true);
		sprintf(name, "min_in_state");
	}

	int vertex_in_layer(int layer);
};



// Vertex that it is in the least number of states
struct RootOrdering : IS_Ordering {

	RootOrdering(IndepSetInst *_inst) : IS_Ordering(_inst, RootOrder) {
		sprintf(name, "root_ordering");
	}

	int vertex_in_layer(BDD* bdd, int layer) { exit(0); }
};



// Vertex that it is in the least number of states
struct LexOrdering : IS_Ordering {

	LexOrdering(IndepSetInst *_inst) : IS_Ordering(_inst, LexOrder) {
		sprintf(name, "lex_ordering");
	}

	int vertex_in_layer(int layer) { exit(0); }
};




// Maximal Path Decomposition
struct MaximalPathDecomp : IS_Ordering {

	vector<int> v_in_layer;   // vertex at each layer

	MaximalPathDecomp(IndepSetInst *_inst) : IS_Ordering(_inst, MaximalPath) {
		sprintf(name, "maxpath");
		construct_ordering();
	}

	int vertex_in_layer(BDD* bdd, int layer) {
		assert( layer >= 0 && layer < inst->graph->n_vertices);
		return v_in_layer[layer];
	}

private:
	void construct_ordering();
};



// Minimum degree ordering
struct MinDegreeOrdering : IS_Ordering {

  vector<int> v_in_layer;   // vertex at each layer

  MinDegreeOrdering(IndepSetInst *_inst) : IS_Ordering(_inst, MinDegree) {
    sprintf(name, "mindegree");
    construct_ordering();
  }

  int vertex_in_layer(BDD* bdd, int layer) {
    assert( layer >= 0 && layer < inst->graph->n_vertices);
    return v_in_layer[layer];
  }
  
private:
  void construct_ordering();
};


// Random ordering
struct RandOrdering : IS_Ordering {

	vector<int> v_in_layer;   // vertex at each layer
	RandOrdering(IndepSetInst *_inst) : IS_Ordering(_inst, Random) {
		sprintf(name, "random");
		v_in_layer.resize(inst->graph->n_vertices);
		construct_ordering();
	}

	int vertex_in_layer(BDD* bdd, int layer) {
		assert( layer >= 0 && layer < inst->graph->n_vertices);
		return v_in_layer[layer];
	}

private:
	void construct_ordering();
};

// FixedOrdering: read from a file
struct FixedOrdering : IS_Ordering {

	vector<int> v_in_layer;   // vertex at each layer
	FixedOrdering(IndepSetInst *_inst, int* sol, int act_size) : IS_Ordering(_inst, RootOrder) {

		sprintf(name, "fixed");
		v_in_layer.resize((act_size+inst->graph->n_vertices)*10);
		set_ordering(sol, act_size);
	}

	int vertex_in_layer(int layer) {
		return v_in_layer[layer];
	}


public:
	void set_ordering(int* sol, int act_size);
};




// Cut vertex decomposition
struct CutVertexDecompositionGeneralGraph : IS_Ordering {

	vector<int> v_in_layer;      // vertex at each layer
	bool** original_adj_matrix;

	CutVertexDecompositionGeneralGraph(IndepSetInst *_inst) : IS_Ordering(_inst, CutVertexGen) {
		sprintf(name, "cut-vertex-gen");
		v_in_layer.resize(inst->graph->n_vertices);

		restrict_graph();
		construct_ordering();
		regenerate_graph();
	}

	int vertex_in_layer(BDD* bdd, int layer) {
		assert( layer >= 0 && layer < inst->graph->n_vertices);
		return v_in_layer[layer];
	}

private:
	void        restrict_graph();
	void        regenerate_graph();
	void        construct_ordering();
	void        identify_components(vector< vector<int> > &comps, vector<bool> &is_in_graph);
	vector<int> find_ordering(vector<bool> is_in_graph);
};


// Cut vertex decomposition
struct CutVertexDecomposition : IS_Ordering {

	vector<int> v_in_layer;   // vertex at each layer

	CutVertexDecomposition(IndepSetInst *_inst) : IS_Ordering(_inst, CutVertex) {
		sprintf(name, "cut-vertex");
		v_in_layer.resize(inst->graph->n_vertices);
		construct_ordering();
	}

	int vertex_in_layer(BDD* bdd, int layer) {
		assert( layer >= 0 && layer < inst->graph->n_vertices);
		return v_in_layer[layer];
	}

private:
	void        construct_ordering();
	void        identify_components(vector< vector<int> > &comps, vector<bool> &is_in_graph);
	vector<int> find_ordering(vector<bool> is_in_graph);
};






#endif
