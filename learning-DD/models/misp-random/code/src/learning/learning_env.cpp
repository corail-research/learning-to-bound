/* MIT License

[Initial work] Copyright (c) 2018 Dai, Hanjun and Khalil, Elias B and Zhang, Yuyu and Dilkina, Bistra and Song, Le
[Adaptation] Copyright (c) 2018 Quentin Cappart, Emmanuel Goutierre, David Bergman and Louis-Martin Rousseau

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
 */

#include "config.h"
#include "learning_env.h"
#include "graph.h"
#include <cassert>
#include <random>
#include <cstdlib>

#include "indepset_solver.hpp"
#include "orderings.hpp"
#include "merge.hpp"
#include "intset.hpp"
#include "instance.hpp"
#include "stats.hpp"


int bdd_max_width = 10000; // -1 if exact
char reward_type = 'W';
char bdd_type = 'U';
double r_scaling = 1;

LearningEnv::LearningEnv() : IEnv() {

}

void LearningEnv::s0(std::shared_ptr<Graph2> _g, bool isTrain) {
    graph = _g;
    covered_set.clear();
    action_list.clear();
    state_seq.clear();
    act_seq.clear();
    reward_seq.clear();
    sum_rewards.clear();
    width = 0;
    bound = 0;

    inst = new IndepSetInst2;
    inst->build_complete_instance(graph->adj_list);

    solver = new IndepSetSolver(inst, bdd_max_width);
    solver->ordering = new OnlineOrdering(inst);
    solver->merger =  new MinLongestPath(inst, bdd_max_width);

    IntSet starting_state;
    starting_state.resize(0, inst->graph->n_vertices-1, true);
    solver->initialize(starting_state, 0);

}

double LearningEnv::step(int a) {
//cout << "step " << a << endl;
    assert(graph);
    assert(covered_set.count(a) == 0);

    state_seq.push_back(action_list);
    act_seq.push_back(a);

    covered_set.insert(a);
    action_list.push_back(a);

    double old_width = width;
    double old_bound = bound;
    double r_t = 0;

    /*if(bdd_type == 'L')
        solver->generate_next_step_restriction(a);
    else if(bdd_type == 'U')
        solver->generate_next_step_relaxation(a);
    else {
        std::cerr << "unknown bdd_type type"  <<  bdd_type << std::endl;
        exit(0);
    }

    width = solver->final_width;
    bound = solver->get_bound();

    if (reward_type == 'W')
        r_t = getReward(old_width);
    else if (reward_type == 'B' && bdd_type == 'U')
        r_t = getRewardUpperBound(old_bound);
    else if (reward_type == 'B' && bdd_type == 'L')
        r_t = getRewardLowerBound(old_bound);
    else if (reward_type == 'M')
        r_t = getRewardMerge();
    else {
        std::cerr << "Unknown reward type"  <<  cfg::reward_type << std::endl;
        exit(0);
    }*/

    reward_seq.push_back(r_t);
    sum_rewards.push_back(r_t);

    return r_t;
}

int LearningEnv::randomAction()
{
    assert(graph);
    avail_list.clear();

    for (int i = 0; i < graph->num_nodes; ++i) {
        if (covered_set.count(i) == 0) {
            avail_list.push_back(i);
        }
    }

    assert(avail_list.size());
    int idx = rand() % avail_list.size();

    return avail_list[idx];
}

void LearningEnv::Insert(IndepSetInst *inst, int g_id, GSet *GSetTest2)
{
	//insert graph from B&B
	Graph* g = inst->graph;
	
	int num_nodes = g->n_vertices;
	int num_edges = g->n_edges;

	cout << "nb nodes " << num_nodes << endl;
	cout << "nb edges " << num_edges << endl;

	bool *present = new bool[num_nodes];
	memset(present, false, sizeof(bool)*num_nodes);

	int *e_from = new int[num_edges];
	int *e_to = new int[num_edges];
	double *w = new double[num_edges];	


	int c = 0;
	for (int i = 0; i < num_nodes; ++i) {
		for (int j = 0; j < g->adj_list[i].size(); ++j) {
			if (i < g->adj_list[i][j]) {
				e_from[c] = i;
				e_to[c] = g->adj_list[i][j];
				w[c] = 1.0;
				//if (i < 1000) {
					//cout << "e " << e_from[c] << " " << e_to[c] << endl;}
				present[i] = true;
				present[g->adj_list[i][j]] = true;
				c++;
			}
		}
	}


	assert(c == num_edges);
	//cout << "num edges checked" << endl;
	
	for (int i = 0; i < num_nodes; ++i) {
		assert(	present[i]);
	}
	//cout << "vertices checked" << endl;

	for (int i = 0; i < num_edges; ++i) {
		assert(e_from[i] < e_to[i]);
		for (int j = 0; j < i; ++j) {
			assert(e_from[i] != e_from[j] || e_to[i] != e_to[j]);
	}}
	//cout << "edges checked" << endl;

	auto g2 = std::make_shared<Graph2>(num_nodes, num_edges, e_from, e_to, w);
	GSetTest2->InsertGraph2(g_id, g2);
	
}

bool LearningEnv::isTerminal() {

    assert(graph);
    return (int) action_list.size() == graph->num_nodes;
}

double LearningEnv::getReward(int old_width) {
    return -r_scaling * (width - old_width); // increase in width is penalized, decrease are rewarded
}

double LearningEnv::getRewardUpperBound(int old_bound) {
    return -r_scaling * (bound - old_bound); // increase in width is penalized, decrease are rewarded
}

double LearningEnv::getRewardLowerBound(int old_bound) {
    return r_scaling * (bound - old_bound); // increase in width is penalized, decrease are rewarded
}

double LearningEnv::getRewardMerge() {
    return -r_scaling * solver->merger->gap;
}
