#!/bin/bash

DATA_PREFIX="${1:-/data/zmen002/kdtree}"
NODE_SIZE="${2:-1000000000}"
DIMENSION="${3:-2}"

./ae_run_incre_update.sh "${DATA_PREFIX}" "${NODE_SIZE}" "${DIMENSION}"
./ae_run_range_query_with_log.sh "${DATA_PREFIX}" "${NODE_SIZE}" "${DIMENSION}"
./ae_run_real_world.sh "${DATA_PREFIX}" "${NODE_SIZE}" "${DIMENSION}"
./ae_scalability.sh "${DATA_PREFIX}" "${NODE_SIZE}" "${DIMENSION}"
