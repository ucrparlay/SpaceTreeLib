#!/bin/bash

DATA_PREFIX="${1:-${DATA_PREFIX:-/data}}"
NODE_SIZE="${2:-${NODE_SIZE:-1000000000}}"

# Configure cmake in build directory
mkdir -p ../build
cd ../build
rm -rf *
cmake -DDEBUG=OFF -DJEMA=ON .. # TODO: make the compiler a parameter as well
make -j
cd ../script_ae

echo "${DATA_PREFIX}"
echo "${NODE_SIZE}"

printf "\n \n =================GENERATE DATA================= \n \n "
./ae_data_generate.sh "${DATA_PREFIX}" "${NODE_SIZE}"

printf "\n \n =================RUN INCRE UPDATE================= \n \n "
mkdir logs
rm -rf logs/*
./ae_run_incre_update.sh "${DATA_PREFIX}" "${NODE_SIZE}" "2"

printf "\n \n =================RUN RANGE QUERY WITH LOG================= \n \n "
./ae_run_range_query_with_log.sh "${DATA_PREFIX}" "${NODE_SIZE}" "2"

printf "\n \n =================RUN REAL WORLD================= \n \n "
./ae_run_real_world.sh "${DATA_PREFIX}" "${NODE_SIZE}" "2"

printf "\n \n =================RUN SCALABILITY================= \n \n "
./ae_scalability.sh "${DATA_PREFIX}" "${NODE_SIZE}" "2"

printf "\n \n =================RUN BATCH UPDATES================= \n \n "
./ae_run_batch_updates.sh "${DATA_PREFIX}" "${NODE_SIZE}" "2"

printf "\n \n =================RUN INCRE UPDATE HIGH DIM================= \n \n "
./ae_run_incre_update_3d.sh "${DATA_PREFIX}" "${NODE_SIZE}" "3"

printf "\n \n =================MERGE INCRE UPDATE================= \n \n "
mkdir -p data
rm -rf data/*
python3 merge_incre.py incre_insert
python3 merge_incre.py incre_delete
python3 merge_range_query.py
python3 merge_real_world.py
python3 merge_scalability.py
python3 merge_batch_updates.py
python3 merge_incre_3d.py incre_insert
python3 merge_incre_3d.py incre_delete

printf "\n \n =================PLOT================= \n \n "
mkdir -p plots
rm -rf plots/*
Rscript plot_fig3_answer_table.R
Rscript plot_fig4_knn.R
Rscript plot_fig5_range_query_scatter.R
Rscript plot_fig7_real_world.R
Rscript plot_fig8_batch_updates.R
Rscript plot_fig9_scalability.R
Rscript plot_fig10_3d_summary.R

printf "\n \n =================CLEAN DATA================= \n \n "
rm -rf "${DATA_PREFIX}/ss_varden_bigint"
rm -rf "${DATA_PREFIX}/uniform_bigint"
rm -rf "${DATA_PREFIX}/ss_varden"
rm -rf "${DATA_PREFIX}/uniform"
rm -rf "${DATA_PREFIX}/geometry"

printf "\n \n =================DONE! HAVE A GOOD DAY!================= \n \n "
