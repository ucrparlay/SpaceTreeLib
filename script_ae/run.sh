#!/bin/bash

DATA_PREFIX="${1:-/data/zmen002/kdtree}"
NODE_SIZE="${2:-1000000000}"
DIMENSION="${3:-2}"

# Configure cmake in build directory
mkdir -p ../build
cd ../build
cmake -DDEBUG=OFF .. # TODO: make the compiler a parameter as well
make -j
cd ../script_ae

printf "\n \n =================GENERATE DATA================= \n \n "
./ae_data_generate.sh "${DATA_PREFIX}" "${NODE_SIZE}" "${DIMENSION}"

printf "\n \n =================RUN INCRE UPDATE================= \n \n "
# ./ae_run_incre_update.sh "${DATA_PREFIX}" "${NODE_SIZE}" "${DIMENSION}"

printf "\n \n =================RUN RANGE QUERY WITH LOG================= \n \n "
# ./ae_run_range_query_with_log.sh "${DATA_PREFIX}" "${NODE_SIZE}" "${DIMENSION}"

printf "\n \n =================RUN REAL WORLD================= \n \n "
# ./ae_run_real_world.sh "${DATA_PREFIX}" "${NODE_SIZE}" "${DIMENSION}"

printf "\n \n =================RUN SCALABILITY================= \n \n "
# ./ae_scalability.sh "${DATA_PREFIX}" "${NODE_SIZE}" "${DIMENSION}"

printf "\n \n =================MERGE INCRE UPDATE================= \n \n "
mkdir -p logs
mkdir -p data
rm -rf logs/*
rm -rf data/*
python3 merge_incre.py incre_insert
python3 merge_incre.py incre_delete
python3 merge_range_query.py
python3 merge_real_world.py
python3 merge_scalability.py

printf "\n \n =================PLOT================= \n \n "
mkdir -p plots
rm -rf plots/*
Rscript plot_fig3_answer_table.R
Rscript plot_fig4_knn.R
Rscript plot_fig5_range_query_scatter.R
Rscript plot_fig7_real_world.R
Rscript plot_fig8_batch_updates.R
Rscript plot_fig9_scalability.R
Rscript plot_fig10_3D_summary.R

printf "\n \n =================DONE! HAVE A GOOD DAY!================= \n \n "
