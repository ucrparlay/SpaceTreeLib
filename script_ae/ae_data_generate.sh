#!/bin/bash
#
set -o xtrace

DATA_PREFIX="${1:-/data/zmen002/kdtree}"
NODE_SIZE="${2:-1000000000}"

make -C ../build/ data_generator data_washer

DIMENSION=2
./../build/data_generator -p ${DATA_PREFIX} -n ${NODE_SIZE} -d ${DIMENSION} -file_num 2 -varden 0
./../build/data_generator -p ${DATA_PREFIX} -n ${NODE_SIZE} -d ${DIMENSION} -file_num 2 -varden 1
./../build/data_washer -coord_type 0 -d ${DIMENSION} -usage 2 -p "${DATA_PREFIX}/uniform_bigint/${NODE_SIZE}_${DIMENSION}/1.in" -output_suffix "_sort_by_0.in"
./../build/data_washer -coord_type 0 -d ${DIMENSION} -usage 2 -p "${DATA_PREFIX}/uniform_bigint/${NODE_SIZE}_${DIMENSION}/2.in" -output_suffix "_sort_by_0.in"

DIMENSION=3
./../build/data_generator -p ${DATA_PREFIX} -n ${NODE_SIZE} -d ${DIMENSION} -file_num 2 -varden 0 -axis_max 1000000
./../build/data_generator -p ${DATA_PREFIX} -n ${NODE_SIZE} -d ${DIMENSION} -file_num 2 -varden 1 -axis_max 1000000
./../build/data_washer -coord_type 0 -d ${DIMENSION} -usage 2 -p "${DATA_PREFIX}/uniform/${NODE_SIZE}_${DIMENSION}/1.in" -output_suffix "_sort_by_0.in"
./../build/data_washer -coord_type 0 -d ${DIMENSION} -usage 2 -p "${DATA_PREFIX}/uniform/${NODE_SIZE}_${DIMENSION}/2.in" -output_suffix "_sort_by_0.in"
