#!/bin/bash

source config_run.sh

file="$1"

if [ "$#" -eq 2 ]; then
    bdd_type="$2"
fi




# Folder of the trained model (inferred from the previous parameters)
result_root=results-local/$g_type-$density/nodes-$min_n-$max_n/seed-$seed/bdd_type-$bdd_type/reward-$reward_type-$bdd_max_width
save_dir=$result_root/ntype-$net_type-embed-$embed_dim-nbp-$max_bp_iter-rh-$reg_hidden-decay-$decay-step-$n_step-batch_size-$batch_size-r_scaling-$r_scaling

# Others
dev_id=0 # gpu card id

python misp_evaluate_dimacs.py \
        -net_type $net_type \
        -n_step $n_step \
        -dev_id $dev_id \
        -decay $decay \
        -test_min_n $test_min_n \
        -test_max_n $test_max_n \
        -test_g_type $test_g_type \
        -test_density $test_density \
        -num_env $num_env \
        -max_iter $max_iter \
        -mem_size $mem_size \
        -learning_rate $learning_rate \
        -max_bp_iter $max_bp_iter \
        -net_type $net_type \
        -max_iter $max_iter \
        -save_dir $save_dir \
        -embed_dim $embed_dim \
        -batch_size $batch_size \
        -reg_hidden $reg_hidden \
        -momentum 0.9 \
        -l2 0.00 \
        -seed $test_seed \
        -w_scale $w_scale \
        -test_min_n $test_min_n \
        -test_max_n $test_max_n \
        -r_scaling $r_scaling \
        -reward_type $reward_type \
        -bdd_type $bdd_type \
        -bdd_max_width $bdd_max_width \
    	-path "$file"
