/*
 * Taken from http://www.andrew.cmu.edu/user/vanhoeve/mdd/ with the notice:
 * "The software can be freely used but comes with no warranty"
 * --------------------------------------------------------
 * Ordering class
 * --------------------------------------------------------
 */

#ifndef ORDERING_HPP2_
#define ORDERING_HPP2_

#include <cassert>
#include <cstdio>
#include "instance.hpp"
#include "bdd2.hpp"
#include "../../../../../../indepset/sources/code/src/indepset_orderings.hpp"
#include <vector>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>


struct IntComparator2 {
	vector<int> &v;
	IntComparator2(vector<int> &_v) : v(_v) { }
	bool operator()(int i, int j) {
		return (v[i] > v[j]);
	}
};


using namespace std;

//enum OrderType {
//  MinState, RandMinState, MaximalPath, Random, CutVertexGen, CutVertex, Fixed,
//  MinDegree
//};

// Class representing a general ordering
struct IS_Ordering2 {

	IndepSetInst2   *inst;
	char           name[256];
	OrderType	   order_type;

	IS_Ordering2(IndepSetInst2* _inst, OrderType _order_type) : inst(_inst), order_type(_order_type) { }

	// returns vertex corresponding to particular layer
	virtual int vertex_in_layer(BDD2* bdd, int layer) = 0;
};


// Vertex that it is in the least number of states
struct MinInState2 : IS_Ordering2 {

	vector<bool> avail_v;   // available vertices

	MinInState2(IndepSetInst2 *_inst) : IS_Ordering2(_inst, MinState) {
		avail_v.resize(inst->graph->n_vertices, true);
		sprintf(name, "min_in_state");
	}

	int vertex_in_layer(BDD2* bdd, int layer);
};


// Vertex that it is in the least number of states - randomized!!
struct RandomizedMinInState22 : IS_Ordering2 {

	vector<bool> avail_v;   // available vertices
	double prob;			// probability

	boost::random::mt19937 gen;

	RandomizedMinInState22(IndepSetInst2 *_inst, double _prob) : IS_Ordering2(_inst, RandMinState), prob(_prob) {
		avail_v.resize(inst->graph->n_vertices, true);
		sprintf(name, "rand_min");

		gen.seed(inst->graph->n_vertices + inst->graph->n_edges);
	}

	int vertex_in_layer(BDD2* bdd, int layer);
};


// Maximal Path Decomposition
struct MaximalPathDecomp2 : IS_Ordering2 {

	vector<int> v_in_layer;   // vertex at each layer

	MaximalPathDecomp2(IndepSetInst2 *_inst) : IS_Ordering2(_inst, MaximalPath) {
		sprintf(name, "maxpath");
		construct_ordering();
	}

	int vertex_in_layer(BDD2* bdd, int layer) {
		assert( layer >= 0 && layer < inst->graph->n_vertices);
		return v_in_layer[layer];
	}

private:
	void construct_ordering();
};



// Minimum degree ordering
struct MinDegreeOrdering2 : IS_Ordering2 {

  vector<int> v_in_layer;   // vertex at each layer

  MinDegreeOrdering2(IndepSetInst2 *_inst) : IS_Ordering2(_inst, MinDegree) {
    sprintf(name, "mindegree");
    construct_ordering();
  }

  int vertex_in_layer(BDD2* bdd, int layer) {
    assert( layer >= 0 && layer < inst->graph->n_vertices);
    return v_in_layer[layer];
  }
  
private:
  void construct_ordering();
};


// Random ordering
struct RandomOrdering : IS_Ordering2 {

	vector<int> v_in_layer;   // vertex at each layer
	RandomOrdering(IndepSetInst2 *_inst) : IS_Ordering2(_inst, Random) {
		sprintf(name, "random");
		v_in_layer.resize(inst->graph->n_vertices);
		construct_ordering();
	}

	int vertex_in_layer(BDD2* bdd, int layer) {
		assert( layer >= 0 && layer < inst->graph->n_vertices);
		return v_in_layer[layer];
	}

private:
	void construct_ordering();
};



// FixedOrdering2: read from a file
struct FixedOrdering2 : IS_Ordering2 {

	vector<int> v_in_layer;   // vertex at each layer
	FixedOrdering2(IndepSetInst2 *_inst, char* filename) : IS_Ordering2(_inst, Fixed) {

		sprintf(name, "fixed");
		v_in_layer.resize(inst->graph->n_vertices);
		read_ordering(filename);
	}

	int vertex_in_layer(BDD2* bdd, int layer) {
		assert( layer >= 0 && layer < inst->graph->n_vertices);
		return v_in_layer[layer];
	}

private:
	void read_ordering(char* filename);
};

// FixedOrdering2: read from a file
struct FixedOrdering2FromSeq : IS_Ordering2 {

	vector<int> v_in_layer;   // vertex at each layer
	FixedOrdering2FromSeq(IndepSetInst2 *_inst, vector<int> action_list) : IS_Ordering2(_inst, Fixed) {

		sprintf(name, "fixed");
		v_in_layer.resize(inst->graph->n_vertices);

		int i = 0;
		for(auto& a : action_list) {
			v_in_layer[i] = inst->node_mapping.at(a);
			i++;
		}
	}

	int vertex_in_layer(BDD2* bdd, int layer) {
		assert( layer >= 0 && layer < inst->graph->n_vertices);
		return v_in_layer[layer];
	}

};

// FixedOrdering2: read from a file
struct OnlineOrdering : IS_Ordering2 {

	vector<int> v_in_layer;   // vertex at each layer
	OnlineOrdering(IndepSetInst2 *_inst) : IS_Ordering2(_inst, Fixed) {
		sprintf(name, "online");
		v_in_layer.resize(inst->graph->n_vertices);

	}

	int vertex_in_layer(BDD2* bdd, int layer) {
		assert( layer >= 0 && layer < inst->graph->n_vertices);
		return v_in_layer[layer];
	}

};

// Cut vertex decomposition
struct CutVertexDecomposition2GeneralGraph2 : IS_Ordering2 {

	vector<int> v_in_layer;      // vertex at each layer
	bool** original_adj_matrix;

	CutVertexDecomposition2GeneralGraph2(IndepSetInst2 *_inst) : IS_Ordering2(_inst, CutVertexGen) {
		sprintf(name, "cut-vertex-gen");
		v_in_layer.resize(inst->graph->n_vertices);

		restrict_graph();
		construct_ordering();
		regenerate_graph();
	}

	int vertex_in_layer(BDD2* bdd, int layer) {
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
struct CutVertexDecomposition2 : IS_Ordering2 {

	vector<int> v_in_layer;   // vertex at each layer

	CutVertexDecomposition2(IndepSetInst2 *_inst) : IS_Ordering2(_inst, CutVertex) {
		sprintf(name, "cut-vertex");
		v_in_layer.resize(inst->graph->n_vertices);
		construct_ordering();
	}

	int vertex_in_layer(BDD2* bdd, int layer) {
		assert( layer >= 0 && layer < inst->graph->n_vertices);
		return v_in_layer[layer];
	}

private:
	void        construct_ordering();
	void        identify_components(vector< vector<int> > &comps, vector<bool> &is_in_graph);
	vector<int> find_ordering(vector<bool> is_in_graph);
};






#endif
