#!/bin/bash

DATA_PREFIX="${1:-/data/zmen002/kdtree}"

./ae_run_incre_update.sh "${DATA_PREFIX}"
./ae_run_range_query_with_log.sh "${DATA_PREFIX}"
./ae_run_real_world.sh "${DATA_PREFIX}"
./ae_scalability.sh "${DATA_PREFIX}"
