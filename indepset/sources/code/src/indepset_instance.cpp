/**
 * -------------------------------------------------
 * Independent Set structure - Implementation
 * -------------------------------------------------
 */

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "indepset_instance.hpp"
#include "util.hpp"


using namespace std;


//
// Read graph in DIMACS format
//
void Graph::read_dimacs(const char* filename) {

	string buffer;
	char   command;

	ifstream input(filename);

	if (!input.is_open()) {
		cerr << "Error: could not open DIMACS graph file " << filename << endl << endl;
		exit(1);
	}

	int read_edges = 0;
	n_edges = -1;

	int source, target;

	while ( read_edges != n_edges && !input.eof() ) {

		input >> command;

		if (command == 'c') {
			// read comment
			getline(input, buffer);

		} else if (command == 'n') {
			// read weight
			input >> source;
			source--;
			input >> weights[source];

		} else if (command == 'p') {
			// read 'edge' or 'col'
			input >> buffer;

			// read number of vertices and edges
			input >> n_vertices;
			input >> n_edges;

			// allocate adjacent matrix
			adj_m = new bool*[n_vertices];
			for (int i = 0; i < n_vertices; i++) {
				adj_m[i] = new bool[n_vertices];
				memset(adj_m[i], false, sizeof(bool)*n_vertices);
			}

			// allocate adjacent list
			adj_list.resize(n_vertices);

			// initialize weights
			weights = new double[n_vertices];
			for (int i = 0; i < n_vertices; ++i) {
				weights[i] = 1.0;
			}

		} else if (command == 'e') {

			if (input.eof()) {
				break;
			}

			// read edge
			input >> source;
			source--;

			input >> target;
			target--;

			set_adj(source, target);
			set_adj(target, source);

			read_edges++;
		}

	}

	input.close();

	int count_edges = 0;
	for (int i = 0; i < n_vertices; i++) {
		for (int j = i+1; j < n_vertices; j++) {
			if (is_adj(i,j)) {
				count_edges++;
			}
		}
	}

	n_edges = count_edges;
}




//
// Export to gml
//
void Graph::export_to_gml(const char* output) {

	ofstream file(output);
	file << "graph [\n";

	for (int i = 0; i < n_vertices; i++) {
		file << "node [\n";
		file << "\tid " << i << "\n";
		file << "\tlabel \"" << i << "\"\n";

		file << "\t graphics [ \n";
		file << "\t\t type \"ellipse\"\n";
		file << "\t\t hasFill 0 \n";
		file << "\t\t ] \n";

		file << "\t]\n" << endl;
	}
	int total_edges = 0;
	for (int i = 0; i < n_vertices; ++i) {
		for (int j = i+1; j < n_vertices; ++j) {
			if ( !is_adj(i, j) )
				continue;
			file << "edge [\n";
			file << "\t source " << i << "\n";
			file << "\t target " << j << "\n";
			file << "\t]\n";
			total_edges++;
		}
	}
	file << "\n]";
	file.close();
	cout << "TOTAL EDGES: " << total_edges << endl;
}


//
// Create an isomorphic graph according to a vertex mapping
// Mapping description: mapping[i] = position where vertex i is in new ordering
//
//
Graph::Graph(Graph* graph, vector<int>& mapping)
: n_vertices(graph->n_vertices), n_edges(graph->n_edges)
{
	// allocate adjacent matrix
	adj_m = new bool*[n_vertices];
	for (int i = 0; i < n_vertices; ++i) {
		adj_m[i] = new bool[n_vertices];
		memset(adj_m[i], false, sizeof(bool)*n_vertices);
	}

	// allocate adjacent list
	adj_list.resize(n_vertices);

	// construct graph according to mapping
	for (int i = 0; i < graph->n_vertices; ++i) {
		for (vector<int>::iterator it = graph->adj_list[i].begin();
				it != graph->adj_list[i].end();
				it++)
		{
			set_adj(mapping[i], mapping[*it]);
		}
	}
}


//
// Print graph
//
void Graph::print() {
	cout << "Graph" << endl;
	for (int v = 0; v < n_vertices; ++v) {
		if (adj_list[v].size() != 0) {
			cout << "\t" << v << " --> ";
			for (vector<int>::iterator it = adj_list[v].begin();
					it != adj_list[v].end();
					++it)
			{
				cout << *it << " ";
			}
			cout << endl;
		}
	}
	cout << endl;
}
