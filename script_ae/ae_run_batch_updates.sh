#!/bin/bash
#
set -o xtrace

# /usr/bin/drop_caches

DATA_PREFIX="${1:-/data/zmen002/kdtree}"
NODE_SIZE="${2:-1000000000}"
DIMENSION="${3:-2}"

tag=128
dim=2
k=10
summary=0
read_file=0
queryType=0
rounds=3

log_path="logs/batch_updates"
mkdir -p "${log_path}"

Tree=(0 1 2 3 5)
make -C ../build/ kd_test p_test baselines

origin_paths=("${DATA_PREFIX}/ss_varden_bigint/${NODE_SIZE}_$((dim))/1.in" "${DATA_PREFIX}/uniform_bigint/${NODE_SIZE}_${dim}/2_sort_by_0.in" "${DATA_PREFIX}/uniform_bigint/${NODE_SIZE}_${dim}/2.in")
insert_paths=("${DATA_PREFIX}/ss_varden_bigint/${NODE_SIZE}_$((dim))/2.in" "${DATA_PREFIX}/uniform_bigint/${NODE_SIZE}_${dim}/1_sort_by_0.in" "${DATA_PREFIX}/uniform_bigint/${NODE_SIZE}_${dim}/1.in")
for tree in "${Tree[@]}"; do
    if [[ ${tree} -eq 0 ]]; then
        solver="kd_test"
        splits=(0)
    elif [[ ${tree} -eq 1 ]]; then
        solver="kd_test"
        splits=(3)
    elif [[ ${tree} -eq 2 ]]; then
        solver="p_test"
        splits=(1 2)
    elif [[ ${tree} -eq 3 ]]; then
        solver="baselines"
        splits=(1 2)
    elif [[ ${tree} -eq 5 ]]; then
        solver="baselines"
        splits=(2)
    else
        solver="baselines"
        splits=(1)
    fi

    for split in "${splits[@]}"; do
        for j in "${!origin_paths[@]}"; do
            file="${origin_paths[${j}]##*/}"
            dest="${log_path}/${tree}_${split}.log"
            : >"${dest}"
            echo ">>>${dest}"
            exe="../build/${solver}"

            input="${origin_paths[${j}]}"
            insert="${insert_paths[${j}]}"
            echo ${input}
            echo ${insert}

            numactl -i all "${exe}" -p "${input}" -I "${insert}" -k ${k} -t ${tag} -d 2 -r ${rounds} -q 0 -i 1 -s 0 -T ${tree} -l ${split} 2>&1 | tee -a "${dest}"
        done
    done
done
