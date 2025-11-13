#!/bin/bash
#
set -o xtrace

DATA_PREFIX="${1:-/data/zmen002/kdtree}"

make -C ../build/ data_generator data_washer

./data_generator -p ${DATA_PREFIX} -n 1000000000 -d 2 -file_num 2 -varden 0
./data_generator -p ${DATA_PREFIX} -n 1000000000 -d 2 -file_num 2 -varden 1
./data_washer -coord_type 0 -d 2 -usage 2 -p "${DATA_PREFIX}/uniform_bigint/1000000000_2/2.in -output_suffix _sort_by_0.in"

# TODO: add real world data
