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


#include "misp_qnet.h"
#include "graph.h"
#include "config.h"

MISPQNet::MISPQNet() : INet()
{
    inputs["node_feat"] = &m_node_feat;
    inputs["edge_feat"] = &m_edge_feat;
    inputs["label"] = &m_y;
    inputs["graph"] = &graph;
    inputs["act_select"] = &m_act_select;
    inputs["rep_global"] = &m_rep_global;
    cfg::node_dim = 2;
    cfg::edge_dim = 4;
}

void MISPQNet::ActNet()
{
    inputs["node_feat"] = &m_node_feat;
    inputs["edge_feat"] = &m_edge_feat;
    inputs["label"] = &m_y;
    inputs["graph"] = &graph;
    inputs["act_select"] = &m_act_select;
    inputs["rep_global"] = &m_rep_global;
    cfg::node_dim = 2;
    cfg::edge_dim = 4;
}

void MISPQNet::BuildNet()
{
    auto graph = add_const< GraphVar >(fg, "graph", true);//"variable graph is already inserted"
    auto action_select = add_const< SpTensorVar<mode, Dtype> >(fg, "act_select", true);
    auto rep_global = add_const< SpTensorVar<mode, Dtype> >(fg, "rep_global", true);

    auto n2nsum_param = af< Node2NodeMsgPass<mode, Dtype> >(fg, {graph});
    auto e2nsum_param = af< Edge2NodeMsgPass<mode, Dtype> >(fg, {graph});
    auto n2esum_param = af< Node2EdgeMsgPass<mode, Dtype> >(fg, {graph});

    auto subgsum_param = af< SubgraphMsgPass<mode, Dtype> >(fg, {graph}, cfg::avg_global);

    cfg::embed_dim = 64;
auto w_n2l = add_diff<DTensorVar>(model, "input-node-to-latent", {cfg::node_dim, cfg::embed_dim});
    auto w_e2l = add_diff<DTensorVar>(model, "input-edge-to-latent", {cfg::edge_dim, cfg::embed_dim});//erreur valgrind
    auto p_node_conv = add_diff< DTensorVar >(model, "linear-node-conv", {cfg::embed_dim, cfg::embed_dim});
    auto trans_node_1 = add_diff< DTensorVar >(model, "trans-node-1", {cfg::embed_dim, cfg::embed_dim});
    auto trans_node_2 = add_diff< DTensorVar >(model, "trans-node-2", {cfg::embed_dim, cfg::embed_dim});

    std::shared_ptr< DTensorVar<mode, Dtype> > h1_weight, h2_weight, last_w;

    if (cfg::reg_hidden > 0)
    {
        h1_weight = add_diff<DTensorVar>(model, "h1_weight", {2 * cfg::embed_dim, cfg::reg_hidden});
        h2_weight = add_diff<DTensorVar>(model, "h2_weight", {cfg::reg_hidden, 1});
        h2_weight->value.SetRandN(0, cfg::w_scale);
        fg.AddParam(h2_weight);
        last_w = h2_weight;
    } else
    {
        h1_weight = add_diff<DTensorVar>(model, "h1_weight", {2 * cfg::embed_dim, 1});
        last_w = h1_weight;
    }

    w_n2l->value.SetRandN(0, cfg::w_scale);
    w_e2l->value.SetRandN(0, cfg::w_scale);
    p_node_conv->value.SetRandN(0, cfg::w_scale);
    trans_node_1->value.SetRandN(0, cfg::w_scale);
    trans_node_2->value.SetRandN(0, cfg::w_scale);
    h1_weight->value.SetRandN(0, cfg::w_scale);
    fg.AddParam(w_n2l);
    fg.AddParam(w_e2l);
    fg.AddParam(p_node_conv);
    fg.AddParam(trans_node_1);
    fg.AddParam(trans_node_2);
    fg.AddParam(h1_weight);

    auto node_input = add_const< DTensorVar<mode, Dtype> >(fg, "node_feat", true);
    auto edge_input = add_const< DTensorVar<mode, Dtype> >(fg, "edge_feat", true);
    auto label = add_const< DTensorVar<mode, Dtype> >(fg, "label", true);

    auto node_init = af<MatMul>(fg, {node_input, w_n2l});
    auto cur_node_embed = af<ReLU>(fg, {node_init});

    auto edge_init = af< MatMul >(fg, {edge_input, w_e2l});

    int lv = 0;
    while (lv < cfg::max_bp_iter)
    {
        lv++;
        auto msg_linear = af< MatMul >(fg, {cur_node_embed, p_node_conv});
        auto n2e_linear = af< MatMul >(fg, {n2esum_param, msg_linear});
        auto edge_rep = af< ElewiseAdd >(fg, {n2e_linear, edge_init});
        edge_rep = af< ReLU >(fg, {edge_rep});
        auto e2n = af<MatMul>(fg, {e2nsum_param, edge_rep});

        auto node_linear = af< MultiMatMul >(fg, {e2n, trans_node_1, cur_node_embed, trans_node_2});
        cur_node_embed = af<ReLU>(fg, {node_linear});
    }

    auto y_potential = af<MatMul>(fg, {subgsum_param, cur_node_embed});

    // Q func given a
    auto action_embed = af<MatMul>(fg, {action_select, cur_node_embed});
    auto embed_s_a = af< ConcatCols >(fg, {action_embed, y_potential});

    auto last_output = embed_s_a;
    if (cfg::reg_hidden > 0)
    {
        auto hidden = af<MatMul>(fg, {embed_s_a, h1_weight});
        last_output = af<ReLU>(fg, {hidden});
    }
    q_pred = af< MatMul >(fg, {last_output, last_w});

    auto diff = af< SquareError >(fg, {q_pred, label});
    loss = af< ReduceMean >(fg, {diff});

    // q func on all a
    auto rep_y = af<MatMul>(fg, {rep_global, y_potential});
    auto embed_s_a_all = af< ConcatCols >(fg, {cur_node_embed, rep_y});

    last_output = embed_s_a_all;
    if (cfg::reg_hidden > 0)
    {
        auto hidden = af<MatMul>(fg, {embed_s_a_all, h1_weight});
        last_output = af<ReLU>(fg, {hidden});
    }

    q_on_all = af< MatMul >(fg, {last_output, last_w});
}

void MISPQNet::SetupGraph2Input(std::vector<int>& idxes,
                              std::vector< std::shared_ptr<Graph2> >& g_list,
                              std::vector< std::vector<int>* >& covered,
                              const int* actions)
{
    int node_cnt = 0, edge_cnt = 0;
    for (size_t i = 0; i < idxes.size(); ++i)
    {
        auto& g = g_list[idxes[i]];
        node_cnt += g->num_nodes;
        edge_cnt += g->num_edges * 2;
    }
    graph.Resize(idxes.size(), node_cnt);
    node_feat.Reshape({(size_t)node_cnt, (size_t)cfg::node_dim});
    node_feat.Fill(1.0);
    edge_feat.Reshape({(size_t)edge_cnt, (size_t)cfg::edge_dim});
    edge_feat.Fill(1.0);

    if (actions)
    {//if (!std::cout.good()) std::cerr << "Oops!\n";
//std::cout << "SIZES " << idxes.size() << " " << (size_t)node_cnt << std::endl;
        act_select.Reshape({idxes.size(), (size_t)node_cnt});
        act_select.ResizeSp(idxes.size(), idxes.size() + 1);
    } else
    {
//std::cout << "SIZES " << idxes.size() << " " << (size_t)cfg::node_dim << " " << (size_t)node_cnt << " " << (size_t)edge_cnt<< std::endl;//1, 2, 39, 30
        rep_global.Reshape({(size_t)node_cnt, idxes.size()});
        rep_global.ResizeSp(node_cnt, node_cnt + 1);
    }
    node_cnt = 0;
    edge_cnt = 0;
    size_t edge_offset = 0;
    for (size_t i = 0; i < idxes.size(); ++i)
    {
        auto& g = g_list[idxes[i]];
        std::set<int> c;
        for (size_t j = 0; j < covered[idxes[i]]->size(); ++j)
        {
            auto& cc = *(covered[idxes[i]]);
            int n_c = cc[j];
            c.insert(n_c);
            node_feat.data->ptr[cfg::node_dim * (node_cnt + n_c)] = 0.0;
        }

        for (int j = 0; j < g->num_nodes; ++j)
        {
            int x = node_cnt + j;
            graph.AddNode(i, x);
            for (auto& p : g->adj_list[j])
            {
		//std::cout << "adding edge " << x << " " << p.first << std::endl;
                graph.AddEdge(edge_cnt, x, node_cnt + p.first);
                auto* edge_ptr = edge_feat.data->ptr + edge_offset;
                edge_ptr[0] = c.count(j);
                edge_ptr[1] = p.second;
                edge_ptr[2] = c.count(p.first) ^ c.count(j);

                edge_offset += cfg::edge_dim;
                edge_cnt++;
            }
		
            if (!actions)
            {
                rep_global.data->row_ptr[node_cnt + j] = node_cnt + j;
                rep_global.data->val[node_cnt + j] = 1.0;
                rep_global.data->col_idx[node_cnt + j] = i;
            }
        }
	//std::cout << edge_cnt << " edges added" << std::endl;	
        if (actions)
        {
            auto act = actions[idxes[i]];
            assert(act >= 0 && act < g->num_nodes);
            act_select.data->row_ptr[i] = i;
            act_select.data->val[i] = 1.0;
            act_select.data->col_idx[i] = node_cnt + act;
        }
        node_cnt += g->num_nodes;
    }
    assert(edge_offset == edge_feat.shape.Count());
    assert(edge_cnt == (int)graph.num_edges);
    assert(node_cnt == (int)graph.num_nodes);
    if (actions)
    {
        act_select.data->row_ptr[idxes.size()] = idxes.size();
        m_act_select.CopyFrom(act_select);
    } else {
        rep_global.data->row_ptr[node_cnt] = node_cnt;
        m_rep_global.CopyFrom(rep_global);
    }

    m_node_feat.CopyFrom(node_feat);
    m_edge_feat.CopyFrom(edge_feat);
}

void MISPQNet::SetupTrain(std::vector<int>& idxes,
                         std::vector< std::shared_ptr<Graph2> >& g_list,
                         std::vector< std::vector<int>* >& covered,
                         std::vector<int>& actions,
                         std::vector<double>& target)
{
    SetupGraph2Input(idxes, g_list, covered, actions.data());

    y.Reshape({idxes.size(), (size_t)1});
    for (size_t i = 0; i < idxes.size(); ++i)
        y.data->ptr[i] = target[idxes[i]];
    m_y.CopyFrom(y);
}

void MISPQNet::SetupPredAll(std::vector<int>& idxes,
                           std::vector< std::shared_ptr<Graph2> >& g_list,
                           std::vector< std::vector<int>* >& covered)
{
    SetupGraph2Input(idxes, g_list, covered, nullptr);
}

MISPQNet::~MISPQNet(){} 
