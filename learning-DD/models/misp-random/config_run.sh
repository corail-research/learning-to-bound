#!/bin/bash


# Seed for the random generation: ensure that the test set remains the same.
sample_name=random
test_seed=42


# Characterics of the training graphs, must be the same as the training
g_type=barabasi_albert # erdos_renyi, barabasi_albert
density=4 # density for ER, attachment parameter for BA
min_n=80 # minimum number of nodes
max_n=85 # maximum number of nodes
seed=42 # Seed used for the training

# Characterics of the tested graphs
test_g_type=barabasi_albert
test_density=4
test_min_n=80
test_max_n=85

# Characterics of the DDs built during the training, must be the same as the training
reward_type=bound
bdd_type=relaxed
bdd_max_width=2


# Parameters used for the learning, must be the same as the training
r_scaling=0.01 # Reward scaling factor
batch_size=32 # max batch size for training/testing
max_bp_iter=4 # max belief propagation iteration
embed_dim=128 # embedding size
net_type=MISPQNet # Network type
decay=1 # discounting factor
reg_hidden=64 # number of hidden layers
learning_rate=0.0001 # learning rate
w_scale=0.01 # init weights with rand normal(0, w_scale)
n_step=1 # number of steps in Q-learning
num_env=10 # number of environments
mem_size=50000 # size of the store for experience replay
max_iter=200000 # number of iterations for the training
