/*
 * ===============================================================
 * IndepSetBDD - Implementation class
 *
 * TODO: memory for branching pool can be fixed, since there is
 *       a maximum number of nodes we will add at a time
 * ===============================================================
 */

#include <algorithm>
#include <iostream>
#include <string>
#include "indepset_bdd.hpp"
#include "config.h"

extern double oracle_init_time;

using namespace std;


//
// Constructor given a DIMACS instance filename
//
// IndepSetBDD(char* filename, int _maxWidth)
// DDX10:
IndepSetBDD::IndepSetBDD(const int _rootWidth,
						 const int _ddWidth,
						 const string & _instanceName,
						 std::string model,
						 const int argc,
						 char** args,
						 int oracle,
						 int rand,
						 int nqval)
	: Inst( new IndepSetInst ),
	  best_lb(0),
	  isExact(false),
	  isLBUpdated(false),
	  initialPosPool(0)
{
	assert( _ddWidth == INF_WIDTH || _ddWidth > 0 );


	//cout << "Reading graph... " << _instanceName.c_str() << endl;
	Inst->read_DIMACS(_instanceName.c_str());
	in_state_counter = new int[Inst->graph->n_vertices];


	// ============================================================
	// Generate initial relaxation to derive variable ordering
	// TODO: save initial UB
	// ============================================================

	cout << "[BDD] Root node computation..." << endl;

	int* sol = (int*) malloc((11 + Inst->graph->n_vertices) * sizeof(int));
	BandBInit(argc, args, model, Inst, sol);


	// generate relaxation		

	active_vertices.resize(Inst->graph->n_vertices);
	for (int i = 0; i < Inst->graph->n_vertices; ++i) {
		active_vertices[i] = i;
	}
	maxWidth = _rootWidth;

	int ub1;

	if(oracle){
		ub1 = generateRelaxationOracle(0, model, argc, args, nqval);}
	else{
		ub1 = generateRelaxationHeuristic(0,rand);}

	upper_bound = ub1;
	cout << "\tUpper bound: " << ub1 << endl;
	cout << "branch node size " << branch_nodes.size() << endl;

	// check if BDD is already exact
	isExact = branch_nodes.empty();
	if (isExact) {
		cout << "is exact from start" << endl;
		return;
	}

	// initialize branch node pool
	for (vector<BranchNode*>::iterator st = branch_nodes.begin(); st != branch_nodes.end(); ++st) {
		(*st)->relax_ub = MIN(ub1, (*st)->longest_path + (int)(*st)->state.size());
	}

	
	// set new variable ordering
	if(0){
	var_ordering = new FixedOrdering(Inst, sol, Inst->graph->n_vertices);
	root_ordering.resize(Inst->graph->n_vertices);
	for( int v = 0; v < Inst->graph->n_vertices; v++ ) {
		root_ordering[var_ordering->vertex_in_layer(v)] = Inst->graph->n_vertices - v - 1;
	}} else {
	var_ordering = new MinInState(Inst);
	}

	// reinitialize width
	maxWidth = _ddWidth;
	free(sol);
	sol = NULL;
}


//
// Generate MDD relaxation. Returns upper bound.
//
// TODO: prune BDD nodes according to best lower bound (using perhaps some estimate?)
//
int IndepSetBDD::generateRelaxationOracle(int initial_lp, std::string model, const int argc, char** args, int nqval) {

	maxWidth = 2;
	int oracle = 1;
	IndepSetInst* old_inst = Inst;
	int* sol;
	int* real_sol;

	vector<int> orig_act_vert(active_vertices);

	sol = (int*) malloc((11 + (Inst->graph->n_vertices*2)) * sizeof(int));
	real_sol = (int*) malloc((active_vertices.size()*2) * sizeof(int));


	n_var = active_vertices.size();
	int it_n = 0;
	int found = 0;
	int found_it = -1;

	//float ordMargin = 3.0/100;
	float ordMargin = 0;

	if(ordMargin != 0){
	for (deque<vector<int>>::iterator it = orderings.begin(); it != orderings.end(); ++it){

		if(n_var > d_orderings[it_n] - n_var*ordMargin && n_var <= d_orderings[it_n]){

			cout << d_orderings[it_n] << " - ";
			for (vector<int>::iterator ord = it->begin(); ord != it->end(); ++ord) {
				cout << (*ord) << " ";
			}
			cout << endl;

			for(int i = 0; i < n_var; i++){
				for (vector<int>::iterator ord = it->begin(); ord != it->end(); ++ord) {
					if(orig_act_vert[i] == (*ord)){
						found = 1;
						break;
					}
					found = 0;
				}
				if(found == 0){
					cout << active_vertices[i] << " not found" << endl;
					break;
				}
			}
			
			if(found == 0){
				cout << it_n << " not matching" << endl;
			}
			else{
				cout << it_n << " matching !" << endl;
				found_it = it_n;
				//break;
			}
		}
		it_n++;
	}
	}

	it_n = 0;
	if(found_it != -1){
		for (vector<int>::iterator ord = orderings[found_it].begin(); ord != orderings[found_it].end(); ++ord) {
			for(int i = 0; i < n_var; i++){
				if(orig_act_vert[i] == (*ord)){
					active_vertices[it_n] = (*ord);
					it_n++;
				}
			}
		}
	}

	if(found_it == -1){
		ask_oracle_count++;
		time_t init_oracle;
		time_t end_init_oracle;
		time(&init_oracle);

		vector<int> alone_vert;

		Inst = BandBCall(old_inst, sol, real_sol, active_vertices, alone_vert, model, argc, args, nqval);

		time(&end_init_oracle);
		oracle_init_time += (double)(end_init_oracle - init_oracle);
	
		vector<int> computed_sol;

		for(unsigned int i = 0; i < active_vertices.size(); i++){

			computed_sol.push_back(real_sol[i]);
			active_vertices[i] = real_sol[i];
		}
		cout << endl;

		if (orderings.size() == maxLenOrd) {
           		orderings.pop_front();
			d_orderings.pop_front();
        	}
		orderings.push_back(computed_sol);
		d_orderings.push_back(active_vertices.size());

	}

	int act_size = active_vertices.size();
	
	var_ordering = new LexOrdering(Inst);

	vertex_in_layer.clear();

	if ( var_ordering->order_type == RootOrder ) {
		
		cout << active_vertices.size() << " active_vertices.size() --" << endl;
		cout << Inst->graph->n_vertices << " inst->graph->n_vertices --" << endl;
		cout << old_inst->graph->n_vertices << " old_inst->graph->n_vertices --" << endl;
		
		root_ordering.resize(active_vertices.size());

		for( int v = 0; v < active_vertices.size(); v++ ) {
			active_vertices[v] = var_ordering->vertex_in_layer(v);
		}
		cout << endl;

		cout << "Order root complete" << endl;

	} else if( var_ordering->order_type == LexOrder ) {
		;
	} else {
		cout << "Order undefined" << endl;
		exit(0);
	}


	State root_state;
	root_state.resize( old_inst->graph->n_vertices, false );
	std::reverse(active_vertices.begin(),active_vertices.end());
	for (vector<int>::const_iterator it = active_vertices.begin(); it != active_vertices.end(); ++it) {
		std::cout << *it << " 	";
		root_state.set(*it, true);
	}
	cout << endl;

	// create initial BDD node
	Node* initial_node = new Node(root_state, initial_lp, true);
	node_pool.clear();
	node_pool[ &(root_state) ] = initial_node;
	root_node = initial_node;

	BDDNodePool::iterator node_it, existing_node_it;
	Node* node;

	// relaxation control variables
	int current_vertex;
	int layer = 0;
	const int num_active_vertices = active_vertices.size();

	while ( layer < num_active_vertices ) {
		// ==================================
		// 1. Vertex and BDD node selection
		// ==================================

		// select next vertex and update active vertex list
		if (var_ordering->order_type == MinState) {
			current_vertex = choose_next_vertex_min_size_next_layer();

		} else if ( var_ordering->order_type == RootOrder ) {
			current_vertex = active_vertices.back();
			active_vertices.pop_back();


		} else if (var_ordering->order_type == LexOrder) {
			current_vertex = active_vertices.back();
			active_vertices.pop_back();

		} else {
			exit(0);
		}

		vertex_in_layer.push_back( current_vertex );
		assert( current_vertex != -1 );

		// Take nodes from the pool that have the current vertex in their state

		nodes_layer.clear();
		node_it = node_pool.begin();

		while (node_it != node_pool.end())	{

			if (node_it->second->state[current_vertex]) {

				// move from the pool to current layer list
				nodes_layer.push_back(node_it->second);//EXAMIN THIS
				node_pool.erase(node_it++);

			} else {
				// Second case: state does not contain vertex; we can maintain it in the pool
				++node_it;
			}
		}

		// PRINT LAYER
		//cout << "Layer " << layer << " - current vertex: " << current_vertex;
		//cout << " - pool size: " << node_pool.size();
		//cout << " - before merge: " << nodes_layer.size();
		//cout << " - total: " << node_pool.size() + nodes_layer.size();
		//cout << endl;


		// ==================================
		// 2. Node merging
		// ==================================
		

		if( maxWidth != INF && (int)nodes_layer.size() > maxWidth ) {
			mergeLayer(layer, nodes_layer, 1);}


		// ==================================
		// 3. Branching
		// ==================================

		Node* branch_node;
		for (vector<Node*>::iterator it = nodes_layer.begin(); it != nodes_layer.end(); ++it) {
			branch_node = (*it);

			// TODO: CHECK THIS!
			if (!branch_node->state[current_vertex]) {
				continue;
			}

			branch_node->state.set(current_vertex, false);

			// --------------- One arc ---------------------

			node = new Node(branch_node->state,
					branch_node->longest_path + old_inst->graph->weights[current_vertex],
					branch_node->exact);


			node->state &= old_inst->adj_mask_compl[current_vertex];

			int estimate_weight = branch_node->longest_path + old_inst->graph->weights[current_vertex];

			estimate_weight += node->state.count();

			if (estimate_weight > best_lb) {

				// Equivalence test: check if node is in list
				existing_node_it = node_pool.find( &(node->state) );
				if (existing_node_it != node_pool.end()) {
					// node already exists in the pool: update node match
					merge(node, existing_node_it->second, 1);
					node = existing_node_it->second;

				} else {
					node_pool[ &(node->state)] = node;
				}
			} else {
				delete node;
			}



			// --------------- Zero arc ---------------------

			// we can use the branch node for the zero arc, since its state was already updated
			node = branch_node;

			// We create a zero arc only if it can improve best lower bound
			// estimate size of the independent set from this vertex. We create it only if it
			// can improve current best upper bound
			estimate_weight = branch_node->longest_path;
//			for (size_t val = node->state.find_first(); val != state_end; val = node->state.find_next(val)) {
//				estimate_weight += inst->graph->weights[val];
//			}
			estimate_weight += node->state.count();

			if (estimate_weight > best_lb) {
				existing_node_it = node_pool.find( &(node->state) );
				if( existing_node_it != node_pool.end() ) {
					// node already exists in the pool: update node match
					merge(node, existing_node_it->second, 1);
					node = existing_node_it->second;

				} else {
					node_pool[ &(node->state)] = node;

					// update eligibility list
					/*if (var_ordering->order_type == MinState) {
						for (vector<int>::iterator v = active_vertices.begin(); v != active_vertices.end(); ++v) {
							if ((node->state)[*v]) {
								in_state_counter[*v]++;
							}
						}
					}*/
				}
			} else {
				delete node;
			}
		}

		// if the number of nodes in the pool is empty, then we do not need to explore this BDD further
		// since there are no better feasible solutions that can be generated from here
		if (node_pool.empty()) {
cout << "empty" << endl;
			while ((int)branch_nodes.size() > initialPosPool) {
				delete branch_nodes.back();
				branch_nodes.pop_back();
			}

			// reset internal parameters
			isExact = false;

			// no new upper bound was generated
			if ( var_ordering->order_type == LexOrder) {
				if(found_it == -1){
					delete Inst;}
					Inst = old_inst;
				if (oracle){
				free(sol);
				sol = NULL;
				free(real_sol);
				real_sol = NULL;}
				delete var_ordering;
			}
			return INF;
		}
		// go to next layer
		layer++;
	}

	// take info from terminal node
	assert( node_pool.size() > 0 );
	const Node* terminal = node_pool.begin()->second;

	upper_bound = terminal->longest_path;
	isExact = terminal->exact;

	delete terminal;

	// if last node is exact, BDD is exact: update lower bound
	if (isExact) {
		updateLocalLB(upper_bound, true);
	}

	// update branch node pool
	for (vector<BranchNode*>::iterator st = branch_nodes.begin(); st != branch_nodes.end(); ++st) {
		(*st)->relax_ub = MIN(	upper_bound,
								(*st)->longest_path + (*st)->state.size()
								);	
	}

	if ( var_ordering->order_type == LexOrder) {
		if(found_it == -1){
			delete Inst;}
			Inst = old_inst;
		if (oracle){
		free(sol);
		sol = NULL;
		free(real_sol);
		real_sol = NULL;
		delete var_ordering;}
	}
	return upper_bound;
}

int IndepSetBDD::generateRelaxationHeuristic(int initial_lp, int rand) {

	maxWidth = 2;
	var_ordering = new MinInState(Inst);
	vertex_in_layer.clear();
	if (var_ordering->order_type == MinState) {
		for (vector<int>::iterator v = active_vertices.begin(); v != active_vertices.end(); ++v) {
			in_state_counter[*v] = 1;
		}
		if(rand){
		cout << "rand ordering" << endl;
		vector<int> v_in_l;
		v_in_l.resize(Inst->graph->n_vertices);
		for( int i = 0; i < Inst->graph->n_vertices; i++ ) {
			v_in_l[i] = i;
		}
		random_shuffle(v_in_l.begin(), v_in_l.end());
		ComparatorAuxIntVectorDescending comp(v_in_l);
		sort( active_vertices.begin(), active_vertices.end(), comp );}

	} else if ( var_ordering->order_type == RootOrder ) {
		ComparatorAuxIntVectorDescending comp(root_ordering);
		sort( active_vertices.begin(), active_vertices.end(), comp );

	} else if ( var_ordering->order_type == Random ) {
		cout << "rand ordering successful" << endl;

	}else if( var_ordering->order_type == LexOrder ) {
		sort( active_vertices.begin(), active_vertices.end() );

	} else {
		cout << "Order undefined" << endl;
		exit(0);
	}


	// ---------------------------------------------------------------------
	// 2. Relaxation
	// ---------------------------------------------------------------------

	// Initialize pool for root node

	// create root node state
	State root_state;
	root_state.resize( Inst->graph->n_vertices, false );
	for (vector<int>::const_iterator it = active_vertices.begin(); it != active_vertices.end(); ++it) {
		root_state.set(*it, true);
	}

	// create initial BDD node
	Node* initial_node = new Node(root_state, initial_lp, true);
	node_pool.clear();
	node_pool[ &(root_state) ] = initial_node;
	root_node = initial_node;

	BDDNodePool::iterator node_it, existing_node_it;
	Node* node;

	// relaxation control variables
	int current_vertex;
	int layer = 0;
	const int num_active_vertices = active_vertices.size();


	while ( layer < num_active_vertices ) {

		// ==================================
		// 1. Vertex and BDD node selection
		// ==================================

		// select next vertex and update active vertex list
		if (var_ordering->order_type == MinState) {
			if(rand){
			current_vertex = active_vertices.back();
			active_vertices.pop_back();}
			else{current_vertex = choose_next_vertex_min_size_next_layer();}

		} else if ( var_ordering->order_type == RootOrder ) {
			current_vertex = active_vertices.back();
			active_vertices.pop_back();


		} else if (var_ordering->order_type == LexOrder) {
			current_vertex = active_vertices.back();
			active_vertices.pop_back();

		} else {
			exit(0);
		}

		vertex_in_layer.push_back( current_vertex );
		assert( current_vertex != -1 );

		// Take nodes from the pool that have the current vertex in their state

		nodes_layer.clear();
		node_it = node_pool.begin();

		while (node_it != node_pool.end())	{

			if (node_it->second->state[current_vertex]) {

				// First case: state contains current vertex, thus we need to branch on it

				// if var ordering is min-in-state, decrement counter of active vertex list
				if (var_ordering->order_type == MinState) {
					for (vector<int>::iterator v = active_vertices.begin(); v != active_vertices.end(); ++v) {
						if ((*node_it->first)[*v]) {
							in_state_counter[*v]--;
						}
					}
				}

				// move from the pool to current layer list
				nodes_layer.push_back(node_it->second);
				node_pool.erase(node_it++);

			} else {
				// Second case: state does not contain vertex; we can maintain it in the pool
				++node_it;
			}
		}

		// PRINT LAYER
//		cout << "Layer " << layer << " - current vertex: " << current_vertex;
//		cout << " - pool size: " << node_pool.size();
//		cout << " - before merge: " << nodes_layer.size();
//		cout << " - total: " << node_pool.size() + nodes_layer.size();
//		cout << endl;


		// ==================================
		// 2. Node merging
		// ==================================

		if( maxWidth != INF && (int)nodes_layer.size() > maxWidth ) {
			mergeLayer(layer, nodes_layer, 0);
		}


		// ==================================
		// 3. Branching
		// ==================================

		Node* branch_node;
		for (vector<Node*>::iterator it = nodes_layer.begin(); it != nodes_layer.end(); ++it) {

			branch_node = (*it);

			// TODO: CHECK THIS!
			if (!branch_node->state[current_vertex]) {
				delete branch_node;
				continue;
			}


			branch_node->state.set(current_vertex, false);

			// --------------- One arc ---------------------

			node = new Node(branch_node->state,
					branch_node->longest_path + Inst->graph->weights[current_vertex],
					branch_node->exact);

			// remove adjacency vertices
			node->state &= Inst->adj_mask_compl[current_vertex];

			// estimate size of the independent set from this vertex. We create it only if it
			// can improve current best upper bound
			int estimate_weight = branch_node->longest_path + Inst->graph->weights[current_vertex];
//			for (size_t val = node->state.find_first(); val != state_end; val = node->state.find_next(val)) {
//				estimate_weight += inst->graph->weights[val];
//			}
			estimate_weight += node->state.count();

			if (estimate_weight > best_lb) {

				// Equivalence test: check if node is in list
				existing_node_it = node_pool.find( &(node->state) );
				if (existing_node_it != node_pool.end()) {
					// node already exists in the pool: update node match
					merge(node, existing_node_it->second, 0);
					node = existing_node_it->second;

				} else {
					node_pool[ &(node->state)] = node;

					// update eligibility list
					if (var_ordering->order_type == MinState) {
						for (vector<int>::iterator v = active_vertices.begin(); v != active_vertices.end(); ++v) {
							if ((node->state)[*v]) {
								in_state_counter[*v]++;
							}
						}
					}
				}
			} else {
				delete node;
			}


			// --------------- Zero arc ---------------------

			// we can use the branch node for the zero arc, since its state was already updated
			node = branch_node;

			// We create a zero arc only if it can improve best lower bound
			// estimate size of the independent set from this vertex. We create it only if it
			// can improve current best upper bound
			estimate_weight = branch_node->longest_path;
//			for (size_t val = node->state.find_first(); val != state_end; val = node->state.find_next(val)) {
//				estimate_weight += inst->graph->weights[val];
//			}
			estimate_weight += node->state.count();

			if (estimate_weight > best_lb) {
				existing_node_it = node_pool.find( &(node->state) );
				if( existing_node_it != node_pool.end() ) {
					// node already exists in the pool: update node match
					merge(node, existing_node_it->second, 0);
					node = existing_node_it->second;

				} else {
					node_pool[ &(node->state)] = node;

					// update eligibility list
					if (var_ordering->order_type == MinState) {
						for (vector<int>::iterator v = active_vertices.begin(); v != active_vertices.end(); ++v) {
							if ((node->state)[*v]) {
								in_state_counter[*v]++;
							}
						}
					}
				}
			} else {
				delete node;
			}
		}

		// if the number of nodes in the pool is empty, then we do not need to explore this BDD further
		// since there are no better feasible solutions that can be generated from here
		if (node_pool.empty()) {

			while ((int)branch_nodes.size() > initialPosPool) {
				delete branch_nodes.back();
				branch_nodes.pop_back();
			}

//			// clean branching pool
//			for (vector<BranchNode*>::iterator it = branch_nodes.begin(); it != branch_nodes.end(); ++it) {
//				delete (*it);
//			}
//			branch_nodes.clear();

			// reset internal parameters
			isExact = false;

			// no new upper bound was generated
			return INF;
		}

		// go to next layer
		layer++;
	}

	// take info from terminal node
	assert( node_pool.size() > 0 );
	const Node* terminal = node_pool.begin()->second;

	upper_bound = terminal->longest_path;
	isExact = terminal->exact;

	delete terminal;
	delete var_ordering;

	// if last node is exact, BDD is exact: update lower bound
	if (isExact) {
		updateLocalLB(upper_bound, true);
	}

	// update branch node pool
	for (vector<BranchNode*>::iterator st = branch_nodes.begin(); st != branch_nodes.end(); ++st) {
		(*st)->relax_ub = MIN(	upper_bound,
								(*st)->longest_path + (*st)->state.size()
								);
	}

	return upper_bound;
}


//
// Merge nodes in a layer to meet maximum allowed width
//
void IndepSetBDD::mergeLayer(int layer, vector<Node*> &nodes_layer, int oa) {

	sort(nodes_layer.begin(), nodes_layer.end(), CompareNodesLongestPath());

	// All nodes will be merged to the one with index "width-1", which will
	// be denoted here by central node
	Node* central_node = nodes_layer[maxWidth-1];

	// If central node is exact, we add it to the branch node pool
	// (since it will be maintained in the BDD, we have to make a copy)
	if( central_node->exact ) {
		addBranchNode(central_node, oa);
		central_node->exact = false;
	}

	// Merge nodes into the central node state (or simply central state)

	State* central_state = &( central_node->state );
	for (vector<Node*>::iterator node = nodes_layer.begin() + maxWidth; node != nodes_layer.end(); ++node) {

		// update state
		(*central_state) |= (*node)->state;

		// if any of the remaining nodes to merge is exact, we have to add it to the branch pool as well.
		if ((*node)->exact) {
			addBranchNode((*node), oa);
		}

		// delete node
		delete (*node);
	}

	// resize node layer vector
	nodes_layer.resize(maxWidth);

	// Check if there are any other nodes which are equivalent to the central node

	Node* node;
	for (int i = 0; i <= maxWidth-2; ++i) {
		node = nodes_layer[i];

		// check if this state already exists in layer nodes
		if( node->state == (*central_state) ) {

			// notice that, if node exists, we do not need to update
			// the costs because the vector is already ordered

			// if existing node is exact, we have to add it to the branch pool
			// (since it will be maintained in the BDD, we have to make a copy)
			if (node->exact) {
				addBranchNode(node, oa);
				node->exact = false;
			}

			// delete the last node (i.e., the central one)
			delete nodes_layer.back();

			// remove it from queue
			nodes_layer.pop_back();
			break;
		}
	}
}


//
// From min-in-state ordering: choose vertex that participates in the least number of states
// and remove it from the active vertex list
//
int IndepSetBDD::choose_next_vertex_min_size_next_layer() {
	int sel_index = 0;
	for (size_t i = 1; i < active_vertices.size(); ++i) {

		if (in_state_counter[active_vertices[i]] != 0 && in_state_counter[active_vertices[i]] < in_state_counter[active_vertices[sel_index]]) {
			sel_index = i;

		} else if (in_state_counter[active_vertices[i]] != 0 && (in_state_counter[active_vertices[i]] == in_state_counter[active_vertices[sel_index]]) && active_vertices[i] < active_vertices[sel_index]) {
			// lexicographic tie breaking
			sel_index = i;
		}
	}

	int sel_vertex = active_vertices[sel_index];
	active_vertices[sel_index] = active_vertices.back();
	active_vertices.pop_back();
	return sel_vertex;
}


//
// Generate MDD restriction. Returns lower bound.
//
int IndepSetBDD::generateRestriction(const int initial_lp) {

	maxWidth = 2;
	var_ordering = new MinInState(Inst);

	// ---------------------------------------
	// 1. Initialization
	// ---------------------------------------

	// width is the number of vertices
  //maxWidth = active_vertices.size();
  //maxWidth = MAX(10, 2000.0 / active_vertices.size());
	// initialize orderings
	if( var_ordering->order_type == MinState ) {
		for( vector<int>::iterator v = active_vertices.begin();
				v != active_vertices.end(); ++v )
		{
			in_state_counter[*v] = 1;
		}

	} else if ( var_ordering->order_type == RootOrder ) {
		ComparatorAuxIntVectorDescending comp(root_ordering);
		sort( active_vertices.begin(), active_vertices.end(), comp );

	} else {
		cout << "Order undefined" << endl;
		exit(0);
	}


	// ---------------------------------------
	// 2. Restriction
	// ---------------------------------------

	// create root node state
	State root_state;
	root_state.resize( Inst->graph->n_vertices, false );
	for (vector<int>::const_iterator it = active_vertices.begin(); it != active_vertices.end(); ++it) {
		root_state.set(*it, true);
	}

	// create initial BDD node
	Node* initial_node = new Node(root_state, initial_lp, true);
	node_pool.clear();
	node_pool[ &(root_state) ] = initial_node;

	BDDNodePool::iterator node_it, existing_node_it;
	Node* node;

	// restriction control variables
	int current_vertex;
	int layer = 0;
	const int num_active_vertices = active_vertices.size();

	while ( layer < num_active_vertices ) {
		// select next vertex
		if( var_ordering->order_type == MinState ) {
			current_vertex = choose_next_vertex_min_size_next_layer();

		} else if ( var_ordering->order_type == RootOrder ) {
			current_vertex = active_vertices.back();
			active_vertices.pop_back();

		} else {
			exit(0);
		}

		// Take nodes from the pool that have the current vertex in their state
		nodes_layer.clear();
		node_it = node_pool.begin();
		while( node_it != node_pool.end() )	{

			if( node_it->second->state[current_vertex] ) {

				// if var ordering is min-in-state, decrement active list
				if( var_ordering->order_type == MinState ) {
					for( vector<int>::iterator v = active_vertices.begin();
							v != active_vertices.end();
							++v )
					{
						if( (*node_it->first)[*v] ) {
							in_state_counter[*v]--;
						}
					}
				}

				// move from the pool to current layer list
				nodes_layer.push_back(node_it->second);
				node_pool.erase(node_it++);

			} else {
				++node_it;
			}
		}

		// PRINT LAYER
		//cout << "Layer " << layer << " - current vertex: " << current_vertex;
		//cout << " - pool size: " << node_pool.size();
		//cout << " - before merge: " << nodes_layer.size();
		//cout << " - total: " << node_pool.size() + nodes_layer.size();
		//cout << endl;

		// 2. Restriction
		if( maxWidth != INF && (int)nodes_layer.size() > maxWidth ) {
			restrictLayer(layer, nodes_layer);
		}
		//cout << "restricted" << endl;



		// 3. Branching
		Node* branch_node;
		for( vector<Node*>::iterator it = nodes_layer.begin(); it != nodes_layer.end(); ++it ) {

			branch_node = (*it);
			branch_node->state.set(current_vertex, false);

			// --------------- One arc ---------------------

			node = new Node(branch_node->state,
					branch_node->longest_path + Inst->graph->weights[current_vertex],
					branch_node->exact);


			// remove adjacency vertices
			node->state &= Inst->adj_mask_compl[current_vertex];


			existing_node_it = node_pool.find( &(node->state) );

			if( existing_node_it != node_pool.end() ) {

				// node already exists in the pool: update node match
				existing_node_it->second->longest_path = MAX(node->longest_path,
						existing_node_it->second->longest_path);
				delete node;

			} else {
				node_pool[ &(node->state)] = node;

				// update eligibility list
				if( var_ordering->order_type == MinState ) {
					for( vector<int>::iterator v = active_vertices.begin();
							v != active_vertices.end();
							++v )
					{
						if( (node->state)[*v] ) {
							in_state_counter[*v]++;
						}
					}
				}
			}



			// --------------- Zero arc ---------------------

			node = branch_node;  // we can reuse the branch node!

			existing_node_it = node_pool.find( &(node->state) );
			if( existing_node_it != node_pool.end() ) {

				// node already exists in the pool: update node match
				existing_node_it->second->longest_path = MAX(node->longest_path,
						existing_node_it->second->longest_path);

				delete node;

			} else {
				node_pool[ &(node->state)] = node;

				// update eligibility list
				if( var_ordering->order_type == MinState ) {
					for( vector<int>::iterator v = active_vertices.begin();
							v != active_vertices.end();
							++v )
					{
						if( (node->state)[*v] ) {
							in_state_counter[*v]++;
						}
					}
				}
			}

		}

		// go to next layer
		layer++;
	}


	// take terminal node
	Node* terminal = node_pool.begin()->second;
	long int lb = terminal->longest_path;
	delete var_ordering;
	delete terminal;
	return lb;
}

//
// Reduce the size of the layer
//
void IndepSetBDD::restrictLayer(int layer, vector<Node*> &nodes_layer) {
	sort(nodes_layer.begin(), nodes_layer.end(), CompareNodesLongestPath());
	for( vector<Node*>::iterator node = nodes_layer.begin()+maxWidth;
			node != nodes_layer.end(); ++node)
	{
		delete (*node);
	}
	nodes_layer.resize(maxWidth);
}



