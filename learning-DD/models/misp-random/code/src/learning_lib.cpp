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
#include "learning_lib.h"
#include "graph.h"
#include "nn_api.h"
#include "misp_qnet.h"
#include "nstep_replay_mem.h"
#include "simulator.h"
#include "learning_env.h"
#include <random>
#include <algorithm>
#include <string>
#include <cstdlib>
#include <signal.h>
#include <mkl.h>
#include <ctime>

clock_t init_pred;
clock_t end_pred;


using namespace gnn;
#define inf 2147483647/2

void intHandler(int dummy) {
    exit(0);
}

int LoadModel(std::string filename) {
    ASSERT(net, "please init the lib before use");
    net->model.Load(filename);
    return 0;
}

int SaveModel(const char* filename) {
    ASSERT(net, "please init the lib before use");
    net->model.Save(filename);
    return 0;
}


int use_cache = 0;
std::map<std::vector<int>, std::vector<int>> predict_cache;

std::vector< std::vector<double>* > list_pred;
LearningEnv* test_env;

extern "C"
int Init(const int argc, char** argv) {
    signal(SIGINT, intHandler);


	cfg::LoadParams(argc, argv);
    GpuHandle::Init(cfg::dev_id, 1);


    if (!strcmp(cfg::reward_type, "width"))
        reward_type = 'W';

    else if (!strcmp(cfg::reward_type, "bound"))
        reward_type = 'B';

    else if (!strcmp(cfg::reward_type, "merge"))
        reward_type = 'M';
    else {
        std::cerr << "unknown reward type"  <<  cfg::reward_type << std::endl;
        exit(0);
    }

    if (!strcmp(cfg::bdd_type, "relaxed"))
        bdd_type = 'U';

    else if (!strcmp(cfg::bdd_type, "restricted"))
        bdd_type = 'L';
    else {
        std::cerr << "unknown bdd type"  <<  cfg::bdd_type << std::endl;
        exit(0);
    }


    bdd_max_width = cfg::bdd_max_width;
    r_scaling = cfg::r_scaling;

    if (!strcmp(cfg::net_type, "MISPQNet")){
        net = new MISPQNet();
	int life = 42;}
    else {
        std::cerr << "unknown net type: " <<  cfg::net_type << std::endl;
        exit(0);
    }
    net->BuildNet();

	cfg::mem_size = 50000;

    NStepReplayMem::Init(cfg::mem_size);

    Simulator::Init(cfg::num_env);
    for (int i = 0; i < cfg::num_env; ++i)
        Simulator::env_list[i] = new LearningEnv();
    test_env = new LearningEnv();

    list_pred.resize(cfg::batch_size);
    for (int i = 0; i < cfg::batch_size; ++i)
        list_pred[i] = new std::vector<double>(2010);//(cfg::max_n + 10);
    return 0;
}


int Init2(const int argc, char** argv, int gid, std::string mod) {

    list_pred.resize(cfg::batch_size);
    for (int i = 0; i < cfg::batch_size; ++i){
	delete list_pred[i];
        list_pred[i] = new std::vector<double>(2010);}//(cfg::max_n + 10);
    return 0;
}


int UpdateSnapshot() {
    net->old_model.DeepCopyFrom(net->model);
    return 0;
}

int InsertGraph2(bool isTest, const int g_id, const int num_nodes, const int num_edges, const int* edges_from, const int* edges_to, const double* weights) {


    auto g = std::make_shared<Graph2>(num_nodes, num_edges, edges_from, edges_to, weights);


    if (isTest)
        GSetTest.InsertGraph2(g_id, g);
    else
        GSetTrain.InsertGraph2(g_id, g);
    return 0;
}

int ClearTrainGraph2s() {
    GSetTrain.graph_pool.clear();
    return 0;
}

int PlayGame(const int n_traj, const double eps) {
    Simulator::run_simulator(n_traj, eps);
    return 0;
}

ReplaySample sample;
std::vector<double> list_target;
double Fit(const double lr) {
    NStepReplayMem::Sampling(cfg::batch_size, sample);
    bool ness = false;
    for (int i = 0; i < cfg::batch_size; ++i)
        if (!sample.list_term[i]) {
            ness = true;
            break;
        }
    if (ness)
        PredictWithSnapshot(sample.g_list, sample.list_s_primes, list_pred);

    list_target.resize(cfg::batch_size);
    for (int i = 0; i < cfg::batch_size; ++i) {
        double q_rhs = 0;
        if (!sample.list_term[i])
            q_rhs = cfg::decay * max(sample.g_list[i]->num_nodes, list_pred[i]->data());
        q_rhs += sample.list_rt[i];
        list_target[i] = q_rhs;
    }

    return Fit(lr, sample.g_list, sample.list_st, sample.list_at, list_target);
}

template<typename K, typename V>
void print_map(std::map<K,V> &m)
{
	for (auto & pair: m) {
		std::cout << "key ";
		for(auto it = pair.first.begin(); it != pair.first.end(); ++it) {
			cout << *it << " ";}
		std::cout << std::endl;
		std::cout << "order ";
		for(auto it = pair.second.begin(); it != pair.second.end(); ++it) {
			cout << *it << " ";}
		std::cout << std::endl;
		cout << "----------" << endl;
	}
}

void BandBInit(const int argc, char** args, std::string mod, IndepSetInst *inst, int* sol){
int gid = 1;

	Init(argc, args);

	LoadModel(mod);
	std::cout << "model " << mod << std::endl;
}

int gid_incr = 2;
IndepSetInst* new_inst = nullptr;

IndepSetInst* BandBCall(IndepSetInst *inst, int *sol, int *real_sol, vector<int> active_vertices, vector<int> &alone_vert, std::string mod, const int argc, char** args, int nqval){


	NStepReplayMem::Clear();
	if(gid_incr > 2){delete new_inst;}
	IndepSetInst* new_inst = new IndepSetInst();
	vector<int> real_active_vertices; 
	
	new_inst->read_inst(inst, active_vertices, real_active_vertices, alone_vert);
	test_env->Insert(new_inst, gid_incr, &GSetTest);	

	//access graph
	std::cout << "active vert size "<< active_vertices.size() << endl;
	std::cout << "-- graph printing"  << std::endl;
	std::cout << "gid_incr " << gid_incr << std::endl;
	std::cout << GSetTest.Get(gid_incr)->num_nodes << endl;
	std::cout << GSetTest.Get(gid_incr)->num_edges << endl;
	std::cout << "graph printing --"  << std::endl;

	vector<int> ord;
	vector<int> key = active_vertices;
	std::sort(key.begin(), key.end());
	if (use_cache && predict_cache.find(key) != predict_cache.end()){
	   	ord = predict_cache[key];
		for(unsigned int i = 0; i < active_vertices.size(); i++){
			real_sol[i] = ord[i];
		}
	}
	else{

	//time
	time_t init_oracle;
	time_t end_init_oracle;
	time(&init_oracle);

	//asking ordering
	GetSol(gid_incr, sol, &GSetTest, nqval);
	GSetTest.graph_pool.clear();

	time(&end_init_oracle);
	//oracle_asking_time += (double)(end_init_oracle - init_oracle);

	new_inst->correct_labels(inst, active_vertices, real_active_vertices, alone_vert, sol, real_sol);
	



    if(use_cache){
	vector<int> pred_ordering;
	pred_ordering.resize(active_vertices.size());
	for(unsigned int i = 0; i < active_vertices.size(); i++){
		pred_ordering[i] = real_sol[i];
	}
      predict_cache[key] = pred_ordering;
    }
	
}
if(gid_incr >= 2){gid_incr++;}
	
	return new_inst;	
}



double GetResult(const int gid, int* sol) {
    std::vector< std::shared_ptr<Graph2> > g_list(1);
    std::vector< std::vector<int>* > states(1);

    test_env->s0(GSetTest.Get(gid),false); 
    states[0] = &(test_env->action_list);
    g_list[0] = test_env->graph;

    double v = 0;
    int new_action;
    while (!test_env->isTerminal())
    {
        Predict(g_list, states, list_pred);
        auto& scores = *(list_pred[0]);
        new_action = arg_max(test_env->graph->num_nodes, scores.data());
        v += test_env->step(new_action) / r_scaling;
    }

    sol[0] = test_env->width;
    sol[1] = test_env->bound;
    return v;
}





double GetSolOS(const int gid, int* sol, GSet* GSetTest2) {

	list_pred.resize(cfg::batch_size);
	for (int i = 0; i < cfg::batch_size; ++i){
		delete list_pred[i];
		list_pred[i] = new std::vector<double>(2010);}//(cfg::max_n + 10);


    std::vector< std::shared_ptr<Graph2> > g_list(1);
    std::vector< std::vector<int>* > states(1);

    test_env->s0(GSetTest2->Get(gid),false);

    states[0] = &(test_env->action_list);//PRINT
    g_list[0] = test_env->graph;

    double v = 0;
    int new_action;

	mkl_free_buffers();


//Time Predict
	init_pred = clock();

        Predict(g_list, states, list_pred);

	end_pred = clock();
	double time_pred = ((double)(end_pred - init_pred))/CLOCKS_PER_SEC;

        auto& scores = *(list_pred[0]);

	for (int i = 0; i < test_env->graph->num_nodes; ++i) {
		new_action = arg_max(test_env->graph->num_nodes, scores.data());
		scores[new_action] = -inf;
		v += test_env->step(new_action) / r_scaling;
	    }

    sol[0] = test_env->graph->num_nodes;
    sol[1] = test_env->width;
    sol[2] = test_env->bound;

    for (int i = 0; i < test_env->graph->num_nodes; ++i) {
        sol[i + 3] = test_env->action_list[i];
    }
	mkl_free_buffers();

    return -v;
}





double GetSol(const int gid, int* sol, GSet* GSetTest2, int nqval) {

	list_pred.resize(cfg::batch_size);
	for (int i = 0; i < cfg::batch_size; ++i){
		delete list_pred[i];
		list_pred[i] = new std::vector<double>(2010);}//(cfg::max_n + 10);

    std::vector< std::shared_ptr<Graph2> > g_list(1);
    std::vector< std::vector<int>* > states(1);

    test_env->s0(GSetTest2->Get(gid),false);

    states[0] = &(test_env->action_list);//PRINT
    g_list[0] = test_env->graph;


    double v = 0;
    int new_action;
    int deltas_cond = 0;
    vector<double> old_scores;
    while (!test_env->isTerminal())
    {

	mkl_free_buffers();



//Time Predict
	init_pred = clock();

        Predict(g_list, states, list_pred);

	end_pred = clock();
	double time_pred = ((double)(end_pred - init_pred))/CLOCKS_PER_SEC;

        auto& scores = *(list_pred[0]);

	for(int i = 0; i<nqval; ++i){
		new_action = arg_max(test_env->graph->num_nodes, scores.data());
		//std::cout << "action " << new_action << " value " << scores.data()[new_action] <<std::endl;
		if(scores.data()[new_action] > -10){v += test_env->step(new_action) / r_scaling;}
		scores.data()[new_action] = -100;
	}


	old_scores = *(list_pred[0]);
    }

    sol[0] = test_env->graph->num_nodes;
    sol[1] = test_env->width;
    sol[2] = test_env->bound;

    for (int i = 0; i < test_env->graph->num_nodes; ++i) {
        sol[i + 3] = test_env->action_list[i];
    }
	mkl_free_buffers();

    return -v;
}

int ClearMem() {
    NStepReplayMem::Clear();
	for (int i = 0; i < cfg::batch_size; ++i){
		delete list_pred[i];}
    return 0;
}
