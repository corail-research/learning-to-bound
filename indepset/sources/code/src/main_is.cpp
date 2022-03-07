//============================================================================
// Name        : BDD_IndepSet.cpp
//============================================================================

#include <ctime>
#include <chrono>
#include <string>
#include <deque>
#include <iostream>
#include <vector>
#include "indepset_instance.hpp"
#include "indepset_bdd.hpp"

//originals settings
//#define MAX_NODES_MEMORY        1700000000
//#define MIN_NODES_DFS		10000


//last test
//#define MAX_NODES_MEMORY        170000000
//#define MIN_NODES_DFS		10000000

//#define MAX_NODES_MEMORY        17000000
//#define MIN_NODES_DFS         15000000

#define MAX_NODES_MEMORY        5000000000
#define MIN_NODES_DFS		100000
#define TIME_LIMIT		1800
//#define TIME_LIMIT            18000

using namespace std;

//
// Buffer of branching nodes
// (Global to preserve memory)
//
vector<BranchNode*> branchNodeBuffer;


//
// Store branching nodes in file. Returns maximum upper bound of stored node
// TODO: do more efficiently: currently it is removing all nodes from queue and adding
// them back
//
inline int storeNodesInFile(const int num_file, const int size, BranchNodeQueue& queue) {

	assert(size < queue.size());
	int maxUB = 0;

	char filename[256];
	sprintf(filename, "branches/input_%d.dat", num_file);

	// create buffer
	branchNodeBuffer.clear();
	for (size_t i = 0; i < queue.size(); ++i) {
		branchNodeBuffer.push_back( queue.top() );
		queue.pop();
	}

	// store node in file
	ofstream states(filename);
	states << size << endl;
	for (int i = 0; i < size; ++i) {

		// get node
		BranchNode* branchNode = branchNodeBuffer.back();
		branchNodeBuffer.pop_back();
		maxUB = MAX(maxUB, branchNode->relax_ub);

		// add to file
		states << branchNode->relax_ub;
		states << " " << branchNode->longest_path;
		states << " " << branchNode->state.size();
		for (vector<int>::iterator it = branchNode->state.begin(); it != branchNode->state.end(); ++it) {
			states << " " << *it << " ";
		}
		states << endl;

		// delete from memory
		delete branchNode;
	}
	states.close();

	// add remaining nodes back to queue
	for (vector<BranchNode*>::iterator it = branchNodeBuffer.begin(); it != branchNodeBuffer.end(); ++it) {
		queue.push(*it);
	}
	branchNodeBuffer.clear();

	return maxUB;
}


//
// Read branching nodes from file
//
inline void readNodesInFile(const int num_file, BranchNodeQueue& queue) {

	char filename[256];
	sprintf(filename, "branches/input_%d.dat", num_file);
	ifstream states(filename);

	int size;
	states >> size;
	for (int i = 0; i < size; ++i) {
		BranchNode* branchNode = new BranchNode;
		states >> branchNode->relax_ub;
		states >> branchNode->longest_path;

		int state_size;
		states >> state_size;
		branchNode->state.resize(state_size);
		for (int j = 0; j < state_size; ++j) {
			states >> branchNode->state[j];
		}

		queue.push(branchNode);
	}
	states.close();
	std::remove(filename);
}


// -------------------------------------------------------
// Global variables
// -------------------------------------------------------
int 					global_lb;
int 					global_ub;
long int 				nodes_explored;
clock_t 				init_time;
time_t 					init_time2;
double 					oracle_init_time;
double 					oracle_asking_time;
double 					test_time;
clock_t 				maxTime;
vector<BranchNode*>		        dfsBranchNodes;
unsigned long int                       size_pool;


//
// Perform DFS search on nodes in the queue
//
void performDFS(IndepSetBDD& indepset_bdd, BranchNodeQueue& queue, std::string model, const int argc, char** args) {

	cout << "\n\nswitching to DFS..." << endl << endl;
	clock_t current;
	//exit(0);

	while (size_pool > MIN_NODES_DFS) {

		// update bound
		global_ub = MIN(global_ub, queue.top()->relax_ub);

		// take best node
		dfsBranchNodes.push_back( queue.top() );
		queue.pop();
		size_pool -= dfsBranchNodes.back()->state.size();

		// explore nodes to reduce size of the queue
		while (!dfsBranchNodes.empty()) {

			// take branch node from the top
			BranchNode* branch_node = dfsBranchNodes.back();
			dfsBranchNodes.pop_back();

			// update statistics
			nodes_explored++;

			if (nodes_explored % 100 == 0 ) {
			  //cout << "**Explored = " << nodes_explored;
				cout << "\t\tTo Explore = " << (queue.size() + dfsBranchNodes.size());
				cout << "\t\tLB = " << global_lb;
				cout << "\t\tUB = " << global_ub;
				//cout << "\t\tgap = " << ((double)(global_ub - global_lb) / global_ub)*100.0;
				cout << "\t\tpool = " << size_pool << endl;
				cout << endl;
			}

			// explores node if not pruned due to global lower bound
			if (branch_node->relax_ub > global_lb) {

				// explore node
				int ub = indepset_bdd.generateRelaxation(branch_node, model, argc, args, 0, 0, 1);

				// update global lower bound if BDD is exact
				if (indepset_bdd.isBDDExact()) {
					global_lb = std::max(global_lb, ub);
				} else {
					// Primal heuristic **************************************
					int lb = indepset_bdd.generateRestriction(branch_node);
					indepset_bdd.updateLocalLB(lb, true);
					global_lb = std::max(global_lb, lb);
					// *******************************************************
				}

				// add open nodes to pool
				indepset_bdd.addBranchNodesDFS(dfsBranchNodes);
			}

			// erase branch node
			delete branch_node;

			// check time
			current = clock();
			if (current - init_time >= maxTime) {
				for (vector<BranchNode*>::iterator it = dfsBranchNodes.begin(); it != dfsBranchNodes.end(); ++it) {
					queue.push(*it);
				}
				return;
			}
		}

		// if bounds match, optimum has been found
		if (global_lb >= global_ub) {
			global_ub = global_lb;
			return;
		}
	}


	cout << "\n\nswitching back to priority queue..." << endl << endl;


}


//
// Main function
//
int main(int argc, char* argv[]) {

	// TODO: add parameters
	if (argc < 4) {
		cout << "\nUsage: indepset [graph file] [root BDD width] [BDD width]\n\n";
		cout << "\twhere: \n";
		cout << "\t\tgraph file: graph in DIMACS format\n";
		cout << "\t\troot BDD width: BDD width at the root node\n";
		cout << "\t\tBDD width: BDD width at branching nodes\n";
		cout << endl;
		exit(1);
	}

	// read input
	int root_max_width = atoi(argv[2]);
	int max_width = atoi(argv[3]);

	int oracle = atoi(argv[6]);
	int rand = 0;

	int nb_arg = 8;
	
	// ------------------------------------------------------
	// Initialize and process initial root node
	// ------------------------------------------------------

	// initialize time
	maxTime = TIME_LIMIT * CLOCKS_PER_SEC;
	init_time = clock();
	time(&init_time2);
	time_t current2;
	clock_t current;
	oracle_init_time = 0;
	test_time = 0;

	int nqval = atoi(argv[5]);	

	IndepSetBDD indepset_bdd(root_max_width, max_width, argv[1], argv[4], argc-nb_arg,(argv+nb_arg), oracle, rand, nqval);



	// Initialize global bounds	
	if (indepset_bdd.isBDDExact()) {
	  cout << "\troot relaxation is exact" << endl << endl;
	  global_ub = indepset_bdd.getUB();
	  global_lb = global_ub;

	} else {
	  global_lb = 0;
	  global_ub = indepset_bdd.getUB();
	}

	// approximate size of the pool
	size_pool = 0;

	// number of nodes explored
	nodes_explored = 1;

	// Create branching pool
	BranchNodeQueue branch_node_queue;
	cout << "MAIN branch node size " << indepset_bdd.getSizeBranch() << endl;
	indepset_bdd.addBranchNodesQueue(branch_node_queue, size_pool);

	// ------------------------------------------------------
	// Search
	// TODO: apply primal heuristic
	// ------------------------------------------------------

	// repeat until all search space is explored
	while (!branch_node_queue.empty()) {

		// switch temporally to DFS if necessary
		if (size_pool > MAX_NODES_MEMORY) {
			performDFS(indepset_bdd, branch_node_queue, argv[4], argc-nb_arg,(argv+nb_arg));
		}

		// check optimality conditions
		if (branch_node_queue.empty()) {
			cout << "out" << endl;
			break;
		}

		// take branch node from the top
		BranchNode* branch_node = branch_node_queue.top();
		branch_node_queue.pop();

		// update statistics
		nodes_explored++;
		size_pool -= branch_node->state.size();

		if (nodes_explored % 100 == 0 ) {
		  //cout << "Explored = " << nodes_explored;
		  cout << "\t\tTo Explore = " << branch_node_queue.size();
		  cout << "\t\tLB = " << global_lb;
		  cout << "\tUB = " << global_ub;
		  cout << "\tgap = " << ((double)(global_ub - global_lb) / global_ub)*100.0;
		  cout << "\t\tpool = " << size_pool;
		  cout << endl;
		}

		// explores node if not pruned due to global lower bound
		if (branch_node->relax_ub > global_lb) {

			time_t init_relax;
			time_t end_relax;
			time(&init_relax);
			auto start = std::chrono::high_resolution_clock::now();

			// explore node
			int ub = indepset_bdd.generateRelaxation(branch_node, argv[4], argc-nb_arg,(argv+nb_arg), oracle, rand, nqval);
			time(&end_relax);
			test_time += (double)(end_relax - init_relax);
			auto finish = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> elapsed = finish - start;

			// update global lower bound if BDD is exact
			if (indepset_bdd.isBDDExact()) {
				global_lb = std::max(global_lb, ub);

			} else {
				// Primal heuristic **************************************
				int lb = indepset_bdd.generateRestriction(branch_node);
				indepset_bdd.updateLocalLB(lb, true);
				global_lb = std::max(global_lb, lb);
				// *******************************************************
			}

			// add open nodes to pool
			indepset_bdd.addBranchNodesQueue(branch_node_queue, size_pool);
		}

		// erase branch node
		delete branch_node;

		// update global upper bound
		if (!branch_node_queue.empty()) {
			global_ub = MIN(global_ub, branch_node_queue.top()->relax_ub);
		} else {
			global_ub = global_lb;
		}

		// if bounds match, optimum has been found
		if (global_lb >= global_ub) {cout << "out2, lb " << global_lb << " - ub " << global_ub << endl;
		  global_ub = global_lb;
		  break;
		}

		//if (global_ub - global_lb < 9){oracle = 1;}
		//if (nodes_explored == 3){oracle = 0;}

		//if (indepset_bdd.getSizeVert() < 40){oracle = 0;}
		//if (indepset_bdd.getSizeVert() >= 40){oracle = 1;}
		

		if (indepset_bdd.ask_oracle_count > atoi(argv[7])-1){oracle = 0;}

		current = clock();
		if (current - init_time >= maxTime) {
			break;
		}
	}

	cout << "out size " << branch_node_queue.size() << endl;

	cout << endl;
	cout << "Lower bound = " << global_lb << endl;
	cout << "Upper bound = " << global_ub << endl;

	current = clock();
	time(&current2);
	double totalTime = ((double)(current - init_time))/CLOCKS_PER_SEC;
	double totTime = (double)(current2 - init_time2);

	char statfilename[256];
	//sprintf(statfilename, "stats.txt");
	sprintf(statfilename, argv[8]);

	ofstream statfile(statfilename, ios::app);
	statfile << argv[1];
	statfile << "\t\t" << "oracle";
	statfile << "\t" << argv[4];
	statfile << "\t\t" << nodes_explored;
	statfile << "\t" << global_lb;
	statfile << "\t" << global_ub;
	statfile << "\t" << ((double)(global_ub - global_lb) / (double)global_ub)*100.0;
	statfile << "\t" << totalTime;
	statfile << "\t\t" << indepset_bdd.ask_oracle_count;
	//statfile << "\t\t" << totTime;
	//statfile << "\t\t" << oracle_init_time;
	//statfile << "\t\t" << oracle_asking_time;
	statfile << endl;
	statfile.close();
	return 0;
}



