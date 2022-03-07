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


#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <vector>
#include <random>
#include "graph.h"

int arg_max(int n, const double* scores);
int arg_min(int n, const double* scores);
double max(int n, const double* scores);

class IEnv;
class Simulator
{
public:
    static void Init(int _num_env);

    static void run_simulator(int num_seq, double eps);

    static int make_action(int num_nodes, std::vector<double>& scores);

    static std::vector<IEnv*> env_list;
    static std::vector< std::shared_ptr<Graph2> > g_list;
    static std::vector< std::vector<int>* > covered;
    static std::vector< std::vector<double>* > pred;

    static std::default_random_engine generator;
    static std::uniform_real_distribution<double> distribution;
};

#endif