#!/bin/bash
set -o xtrace

DATA_PREFIX="${1:-/data/zmen002/kdtree}"

# /usr/bin/drop_caches

Node=(1000000000)
Tree=(0 1 2 3 4)
# Tree=(3)
Dims=(2)
# Type=(1)
# paths=("${DATA_PREFIX}/ss_varden_bigint/1000000000_2/1.in" "${DATA_PREFIX}/uniform_bigint/1000000000_2/2_sort_by_0.in" "${DATA_PREFIX}/uniform_bigint/1000000000_2/2.in")

k=10
insNum=0
summary=0
read_file=0
queryType=$((2#0)) # 1110000
round=2
log_path="logs/range_query_log"
mkdir -p "${log_path}"

make -C ../build/ kd_test p_test baselines boost_rtree

dim=2
paths=("${DATA_PREFIX}/ss_varden_bigint/1000000000_$((dim))/1.in" "${DATA_PREFIX}/uniform_bigint/1000000000_${dim}/2_sort_by_0.in" "${DATA_PREFIX}/uniform_bigint/1000000000_${dim}/2.in")
for tree in "${Tree[@]}"; do
    if [[ ${tree} -eq 0 ]]; then
        solver="kd_test"
        splits=(0)
    elif [[ ${tree} -eq 1 ]]; then
        solver="kd_test"
        splits=(3)
    elif [[ ${tree} -eq 2 ]]; then
        solver="p_test"
        splits=(1)
    elif [[ ${tree} -eq 3 ]]; then
        solver="baselines"
        splits=(1 2)
    elif [[ ${tree} -eq 4 ]]; then
        solver="boost_rtree"
        splits=(0)
    elif [[ ${tree} -eq 5 ]]; then
        solver="baselines"
        splits=(2)
    else
        echo "ERROR: Invalid tree type ${tree}"
        exit 1
    fi

    for split in "${splits[@]}"; do
        tag=64
        dest="${log_path}/${dim}_${tree}_${split}.log"
        : >"${dest}"
        echo ">>>${dest}"
        exe="../build/${solver}"

        for path in "${paths[@]}"; do
            numactl -i all ${exe} -p ${path} -r ${round} -k ${k} -i ${read_file} -s ${summary} -t ${tag} -d ${dim} -q ${queryType} -T ${tree} -l ${split} 2>&1 | tee -a "${dest}"
        done
    done
done
