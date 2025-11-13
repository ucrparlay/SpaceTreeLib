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
./ae_run_incre_update.sh "${DATA_PREFIX}" "${NODE_SIZE}" "${DIMENSION}"

printf "\n \n =================RUN RANGE QUERY WITH LOG================= \n \n "
# ./ae_run_range_query_with_log.sh "${DATA_PREFIX}" "${NODE_SIZE}" "${DIMENSION}"

printf "\n \n =================RUN REAL WORLD================= \n \n "
# ./ae_run_real_world.sh "${DATA_PREFIX}" "${NODE_SIZE}" "${DIMENSION}"

printf "\n \n =================RUN SCALABILITY================= \n \n "
# ./ae_scalability.sh "${DATA_PREFIX}" "${NODE_SIZE}" "${DIMENSION}"
