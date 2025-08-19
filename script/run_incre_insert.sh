#!/bin/bash
set -o xtrace

/usr/bin/drop_caches

Node=(1000000000)
Tree=(0 1 2)
Dim=(2)
Type=(0 1)
# paths=("/data/zmen002/kdtree/ss_varden_bigint/1000000000_2/1.in" "/data/zmen002/kdtree/uniform_bigint/1000000000_2/2_sort_by_0.in" "/data/zmen002/kdtree/uniform_bigint/1000000000_2/2.in")
paths=("/data/zmen002/kdtree/ss_varden/1000000000_2/1.in" "/data/zmen002/kdtree/uniform/1000000000_2/2_sort_by_0.in" "/data/zmen002/kdtree/uniform/1000000000_2/2.in")

k=10
insNum=1
summary=0
read_file=0
queryType=$((2#0)) # 1110000
round=3
resFile=""

make -C ../build/ -j

for query_type in "${Type[@]}"; do

    if [[ ${query_type} -eq 0 ]]; then
        tag=$((2#1000)) # 1110000
        dest="incre_insert.log"
    else
        tag=$((2#10000)) # 1110000
        dest="incre_delete.log"
    fi
    : >"${dest}"
    echo ">>>${dest}"

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
        fi
        exe="../build/${solver}"

        for split in "${splits[@]}"; do
            for path in "${paths[@]}"; do
                numactl -i all ${exe} -p ${path} -r ${round} -k ${k} -i ${read_file} -s ${summary} -t ${tag} -d 2 -q ${queryType} -T ${tree} -l ${split} 2>&1 | tee -a "${dest}"
            done
        done
    done
done

current_date_time="$(date "+%d %H:%M:%S")"
echo $current_date_time
