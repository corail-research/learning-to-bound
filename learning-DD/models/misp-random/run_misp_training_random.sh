#!/bin/bash

source config_run.sh

# Folder to save the trained model
result_root=results-local/$g_type-$density/nodes-$min_n-$max_n/seed-$seed/bdd_type-$bdd_type/reward-$reward_type-$bdd_max_width
save_dir=$result_root/ntype-$net_type-embed-$embed_dim-nbp-$max_bp_iter-rh-$reg_hidden-decay-$decay-step-$n_step-batch_size-$batch_size-r_scaling-$r_scaling

# Others
dev_id=0 # gpu card id
plot_training=1 # Boolean value: plot the training curve or not

if [ ! -e $save_dir ];
then
    mkdir -p $save_dir
fi



python misp_training_random.py \
    -net_type $net_type \
    -n_step $n_step \
    -dev_id $dev_id \
    -decay $decay \
    -min_n $min_n \
    -max_n $max_n \
    -num_env $num_env \
    -max_iter $max_iter \
    -mem_size $mem_size \
    -g_type $g_type \
    -density $density \
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
    -w_scale $w_scale \
    -seed $seed \
    -reward_type $reward_type \
    -bdd_type $bdd_type \
    -bdd_max_width $bdd_max_width \
    -r_scaling $r_scaling \
    -plot_training $plot_training \
    2>&1 | tee $save_dir/log-training.txt


