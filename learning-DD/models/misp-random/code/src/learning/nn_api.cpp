/* MIT License

Copyright (c) 2018 Dai, Hanjun and Khalil, Elias B and Zhang, Yuyu and Dilkina, Bistra and Song, Le

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
#include "nn_api.h"
#include "graph.h"
#include "tensor/tensor_all.h"
#include <random>
#include <algorithm>
#include <cstdlib>
#include <signal.h>
#include <memory>

using namespace gnn;

#define inf 2147483647/2

INet* net = nullptr;
//INet* net;//ok, init in learning_lib.cpp
//std::shared_ptr<INet> net(new MISPQNet());

std::vector<int> batch_idxes;

void Predict(std::vector< std::shared_ptr<Graph2> >& g_list, std::vector< std::vector<int>* >& covered, std::vector< std::vector<double>* >& pred)
{
    DTensor<CPU, Dtype> output;
    int n_graphs = g_list.size();//1
    for (int i = 0; i < n_graphs; i += cfg::batch_size)
    {
        int bsize = cfg::batch_size;//1
        if (i + cfg::batch_size > n_graphs)
            bsize = n_graphs - i;
        batch_idxes.resize(bsize);
        for (int j = i; j < i + bsize; ++j)
            batch_idxes[j - i] = j;

//std::cout << "BATCH IDXES " << batch_idxes[0] << std::endl; //0
//std::cout << "N_C " << (covered[batch_idxes[0]])[0] << std::endl; //pour comprendre node_feat.data->ptr[cfg::node_dim * (node_cnt + n_c)] = 0.0 dans SetupGraphInput, misp_qnet.cpp
//std::cout << "N_C " << (*covered[0])[1] << std::endl; //seg fault

        net->SetupPredAll(batch_idxes, g_list, covered);
//std::cerr << "INPUT " << net->inputs["graph"] << std::endl;//CLARIFY THIS
        net->fg.FeedForward({net->q_on_all}, net->inputs, Phase::TEST);        
        auto& raw_output = net->q_on_all->value;
        output.CopyFrom(raw_output); //No output
        
        int pos = 0;
        for (int j = i; j < i + bsize; ++j)
        {
            auto& cur_pred = *(pred[j]);
            auto& g = g_list[j];
//net->pr(); net INet exists
//std::cerr << "pred	num nodes " << g->num_nodes << std::endl; //graph well loaded

//for (int k = 0; k < g->num_nodes; ++k)
//            {std::cerr << "pred	edge 0 " << g->adj_list[0][k].first << std::endl;
//}

            for (int k = 0; k < g->num_nodes; ++k)
            {
		//if(/*pos < 1*/1){
		//	std::cerr << "pred " << pos << " " << output.data->ptr[pos] << std::endl;}
                cur_pred[k] = output.data->ptr[pos];
                pos += 1;
            }            
            auto& cur_covered = *(covered[j]);//you don't take an action twice !
            for (auto& k : cur_covered)
                cur_pred[k] = -inf;
        }

        ASSERT(pos == (int)output.shape.Count(), "idxes not match");
    }   
}

void PredictWithSnapshot(std::vector< std::shared_ptr<Graph2> >& g_list, std::vector< std::vector<int>* >& covered, std::vector< std::vector<double>* >& pred)
{
    net->UseOldModel();
    Predict(g_list, covered, pred);
    net->UseNewModel();
}

double Fit(const double lr, std::vector< std::shared_ptr<Graph2> >& g_list, std::vector< std::vector<int>* >& covered, std::vector<int>& actions, std::vector<double>& target)
{   
    Dtype loss = 0;
    int n_graphs = g_list.size();
    for (int i = 0; i < n_graphs; i += cfg::batch_size)
    {
        int bsize = cfg::batch_size;
        if (i + cfg::batch_size > n_graphs)
            bsize = n_graphs - i;

        batch_idxes.resize(bsize);
        for (int j = i; j < i + bsize; ++j)
            batch_idxes[j - i] = j;

        net->SetupTrain(batch_idxes, g_list, covered, actions, target);
        net->fg.FeedForward({net->loss}, net->inputs, Phase::TRAIN);
        net->fg.BackPropagate({net->loss});
        net->learner->cur_lr = lr;
        net->learner->Update();

        loss += net->loss->AsScalar() * bsize;
    }
    
    return loss / g_list.size();
}
