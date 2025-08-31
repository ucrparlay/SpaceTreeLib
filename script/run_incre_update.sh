#!/bin/bash
set -o xtrace

/usr/bin/drop_caches

Node=(1000000000)
# Tree=(0 1 2)
Tree=(5)
Dims=(3)
Type=(0 1)
# paths=("/data/zmen002/kdtree/ss_varden_bigint/1000000000_2/1.in" "/data/zmen002/kdtree/uniform_bigint/1000000000_2/2_sort_by_0.in" "/data/zmen002/kdtree/uniform_bigint/1000000000_2/2.in")

k=10
insNum=0
summary=0
read_file=0
# queryType=$((2#0)) # 1110000
round=3
log_path="logs/incre_update"
mkdir -p "${log_path}"

make -C ../build/ kd_test p_test baselines

for dim in "${Dims[@]}"; do
    if [[ ${dim} -eq 2 ]]; then
        paths=("/data/zmen002/kdtree/ss_varden_bigint/1000000000_$((dim))/1.in" "/data/zmen002/kdtree/uniform_bigint/1000000000_${dim}/2_sort_by_0.in" "/data/zmen002/kdtree/uniform_bigint/1000000000_${dim}/2.in")
    elif [[ ${dim} -eq 3 ]]; then
        paths=("/data/zmen002/kdtree/ss_varden/1000000000_$((dim))/1.in" "/data/zmen002/kdtree/uniform/1000000000_${dim}/2_sort_by_0.in" "/data/zmen002/kdtree/uniform/1000000000_${dim}/2.in")
    fi
    for query_type in "${Type[@]}"; do
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
            elif [[ ${tree} -eq 5 ]]; then
                solver="baselines"
                splits=(2)
            else
                solver="baselines"
                splits=(0)
            fi

            for split in "${splits[@]}"; do
                if [[ ${query_type} -eq 0 ]]; then
                    tag=$((2#1000)) # 1110000
                    dest="${log_path}/incre_insert_${dim}_${tree}_${split}.log"
                else
                    tag=$((2#10000)) # 1110000
                    dest="${log_path}/incre_delete_${dim}_${tree}_${split}.log"
                fi
                : >"${dest}"
                echo ">>>${dest}"
                exe="../build/${solver}"

                for path in "${paths[@]}"; do
                    numactl -i all ${exe} -p ${path} -r ${round} -k ${k} -i ${read_file} -s ${summary} -t ${tag} -d ${dim} -q ${queryType} -T ${tree} -l ${split} 2>&1 | tee -a "${dest}"
                done
            done
        done
    done
done
