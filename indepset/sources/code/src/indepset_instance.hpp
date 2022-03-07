/*
 * --------------------------------------------------------
 * Instance data structures
 * --------------------------------------------------------
 */

#ifndef INSTANCE_HPP_
#define INSTANCE_HPP_

#include <cstring>
#include <fstream>
#include <iostream>
#include <boost/dynamic_bitset.hpp>

using namespace std;


//
// Graph structure
//
struct Graph {

	int                         n_vertices;         /**< |V| */
	int                         n_edges;            /**< |E| */
	double*						weights;			/**< weight of each vertex */

	bool**                      adj_m;              /**< adjacent matrix */
	vector< vector<int> >       adj_list;           /**< adjacent list */


	/** Set two vertices as adjacents */
	void set_adj(int i, int j);

	/** Check if two vertices are adjancent */
	bool is_adj(int i, int j);

	/** Empty constructor */
	Graph();

	/** Create an isomorphic graph according to a vertex mapping */
	Graph(Graph* graph, vector<int>& mapping);

	/** Read graph from a DIMACS format */
	void read_dimacs(const char* filename);

	/** Export to GML format */
	void export_to_gml(const char* output);

	/** Constructor with number of vertices */
	Graph(int num_vertices);

	/** Add edge */
	void add_edge(int i, int j);

	/** Remove edge */
	void remove_edge(int i, int j);

	/** Return degree of a vertex */
	int degree( int v ) { return adj_list[v].size(); }

	/** Print graph */
	void print();
};



//
// Independent set instance structure
//
struct IndepSetInst {

	Graph*              				graph;             	// independent set graph
	vector< boost::dynamic_bitset<> >	adj_mask_compl;	 	// complement mask of adjacencies


	/** Read DIMACS independent set instance */
	void read_DIMACS(const char* filename);
	void read_inst(IndepSetInst* inst,vector<int> active_vertices, vector<int> &real_active_vertices, vector<int> &alone_vert);
	void correct_labels(IndepSetInst* inst, vector<int> active_vertices, vector<int> &real_active_vertices, vector<int> &alone_v, int* sol, int* real_sol);
};



/*
 * -----------------------------------------------------
 * Inline implementations: Graph
 * -----------------------------------------------------
 */


/**
 * Empty constructor
 */
inline Graph::Graph() : n_vertices(0), n_edges(0), weights(NULL), adj_m(NULL) {
}


/**
 * Constructor with number of vertices
 **/
inline Graph::Graph(int num_vertices)
: n_vertices(num_vertices), n_edges(0), weights(NULL)
{
	adj_m = new bool*[ num_vertices ];
	for (int i = 0; i < num_vertices; ++i) {
		adj_m[i] = new bool[ num_vertices ];
		memset( adj_m[i], false, sizeof(bool) * num_vertices );
	}
	adj_list.resize(num_vertices);
}


/**
 * Check if two vertices are adjacent
 */
inline bool Graph::is_adj(int i, int j) {
	assert(i >= 0);
	assert(j >= 0);
	assert(i < n_vertices);
	assert(j < n_vertices);
	return adj_m[i][j];
}


/**
 * Set two vertices as adjacent
 */
inline void Graph::set_adj(int i, int j) {
	assert(i >= 0);
	assert(j >= 0);
	assert(i < n_vertices);
	assert(j < n_vertices);

	// check if already adjacent
	if (adj_m[i][j])
		return;

	// add to adjacent matrix and list
	adj_m[i][j] = true;
	adj_list[i].push_back(j);
}



/**
 * Add edge
 **/
inline void Graph::add_edge(int i, int j) {
	assert(i >= 0);
	assert(j >= 0);
	assert(i < n_vertices);
	assert(j < n_vertices);

	// check if already adjacent
	if (adj_m[i][j])
		return;

	// add to adjacent matrix and list
	adj_m[i][j] = true;
	adj_m[j][i] = true;
	adj_list[i].push_back(j);
	adj_list[j].push_back(i);

	n_edges++;
}

/**
 * Remove edge
 **/
inline void Graph::remove_edge(int i, int j) {
	assert(i >= 0);
	assert(j >= 0);
	assert(i < n_vertices);
	assert(j < n_vertices);

	// check if already adjacent
	if (!adj_m[i][j])
		return;

	// add to adjacent matrix and list
	adj_m[i][j] = false;
	adj_m[j][i] = false;

	for (int v = 0; v < (int)adj_list[i].size(); ++v) {
		if ( adj_list[i][v] == j ) {
			adj_list[i][v] = adj_list[i].back();
			adj_list[i].pop_back();
			break;
		}
	}

	for (int v = 0; v < (int)adj_list[j].size(); ++v) {
		if ( adj_list[j][v] == i ) {
			adj_list[j][v] = adj_list[j].back();
			adj_list[j].pop_back();
			break;
		}
	}
}


/*
 * -----------------------------------------------------
 * Inline implementations: Independent Set
 * -----------------------------------------------------
 */


//
// Read DIMACS independent set instance with no costs
//
inline void IndepSetInst::read_DIMACS(const char *filename) {

	cout << "\nReading instance " << filename << endl;

	// read graph
	graph = new Graph;
	graph->read_dimacs(filename);

	cout << "\tnumber of vertices: " << graph->n_vertices << endl;
	cout << "\tnumber of edges: " << graph->n_edges << endl;

	// create complement mask of adjacencies
	adj_mask_compl.resize(graph->n_vertices);
	for( int v = 0; v < graph->n_vertices; v++ ) {

		adj_mask_compl[v].resize(graph->n_vertices, true);
		for( int w = 0; w < graph->n_vertices; w++ ) {
			if( graph->is_adj(v,w) ) {
				adj_mask_compl[v].set(w, false);
			}
		}

		// we assume here a vertex is adjacent to itself
		adj_mask_compl[v].set(v, false);

	}
	cout << "\tdone.\n" << endl;
}

inline void IndepSetInst::read_inst(IndepSetInst* inst,vector<int> active_vertices, vector<int> &real_active_vertices, vector<int> &alone_vert) {
//A new instance is created by keeping only the active vertices in the old instance
//The vertices have their number changed for their position in the active_vertices list
//Unconnected vertices should be ignored /!\ Return alone_vert to be treated in FixedOrdering


	//Alone vertices elimination
	bool *loneliness = new bool[active_vertices.size()];
	memset(loneliness, true, sizeof(bool)*active_vertices.size());

	alone_vert.resize(active_vertices.size());
	real_active_vertices = active_vertices;

	for (unsigned int i = 0; i < active_vertices.size(); ++i){
		for (unsigned int j = i+1; j < active_vertices.size(); ++j){
			if(inst->graph->adj_m[active_vertices[i]][active_vertices[j]]){
				loneliness[i] = false;
				loneliness[j] = false;
			}
	}}


	int c = 0;
	for (int i = 0; i < active_vertices.size(); ++i){
		if(loneliness[i]){
			real_active_vertices.erase(real_active_vertices.begin() + i - c);
			alone_vert[c] = active_vertices[i];
			c++;
	}}
	alone_vert.resize(c);




	// read graph
	graph = new Graph;

	graph->n_vertices = real_active_vertices.size();
	graph->n_edges = 0;

	graph->weights = new double[graph->n_vertices];
	for (int i = 0; i < graph->n_vertices; ++i) {
		graph->weights[i] = 1.0;
	}

	graph->adj_m = new bool*[graph->n_vertices];
	for (int i = 0; i < graph->n_vertices; i++) {
		graph->adj_m[i] = new bool[graph->n_vertices];
		memset(graph->adj_m[i], false, sizeof(bool)*graph->n_vertices);
	}

	// allocate adjacent list
	graph->adj_list.resize(graph->n_vertices);

	for (unsigned int i = 0; i < real_active_vertices.size(); ++i){
		for (unsigned int j = i+1; j < real_active_vertices.size(); ++j){
			if(inst->graph->adj_m[real_active_vertices[i]][real_active_vertices[j]]){
				graph->set_adj(j, i);
				graph->set_adj(i, j);
				graph->n_edges++;
			}
	}}

	
	// create complement mask of adjacencies
	adj_mask_compl.resize(graph->n_vertices);
	for( int v = 0; v < graph->n_vertices; v++ ) {

		adj_mask_compl[v].resize(graph->n_vertices, true);
		for( int w = 0; w < graph->n_vertices; w++ ) {
			if( graph->is_adj(v,w) ) {
				adj_mask_compl[v].set(w, false);
			}
		}

		// we assume here a vertex is adjacent to itself
		adj_mask_compl[v].set(v, false);

	}
	delete[] loneliness;
}




inline void IndepSetInst::correct_labels(IndepSetInst* inst, vector<int> active_vertices, vector<int> &real_active_vertices, vector<int> &alone_v, int* sol, int* real_sol) {


	for(unsigned int v = 0; v < real_active_vertices.size(); v++ ) {
		real_sol[v + alone_v.size()] = real_active_vertices[sol[v + 3]];
	}
	for(unsigned int v = 0; v < alone_v.size(); v++ ) {
		real_sol[v] = alone_v[v];
	}
}


#endif /* INSTANCE_HPP_ */
