#!/bin/bash
set -o xtrace

# Solvers=("test" "zdtree" "cgal")
DataPath="/data/legacy/data3/zmen002/kdtree/geometry/round_small_coord"
declare -A file2Dims
file2Dims["Cosmo50_round_no_dup"]="3"
file2Dims["GeoLifeNoScale_round_no_dup"]="3"
# file2Dims["osm_round_no_dup"]="2"

tag=32
k=10
onecore=0
readFile=0
summary=0
queryType=$((2#0)) # 1110000
# queryType=1 # 1110000
# Tree=(0 1 2 3 4)
Tree=(5)

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
    elif [[ ${tree} -eq 4 ]]; then
        solver="boost_rtree"
        splits=(0)
    elif [[ ${tree} -eq 5 ]]; then
        solver="baselines"
        splits=(2)
    else
        solver="baselines"
        splits=(0)
    fi

    for split in "${splits[@]}"; do
        log_path="logs/real_world"
        mkdir -p ${log_path}
        dest="${log_path}/${tree}_${split}.log"
        : >${dest}
        echo ">>>${dest}"
        exe="../build/${solver}"

        for filename in "${!file2Dims[@]}"; do
            echo ${filename}

            numactl -i all ${exe} -p "${DataPath}/${filename}.in" -k ${k} -t ${tag} -d ${file2Dims[${filename}]} -r 3 -q ${queryType} -i ${readFile} -s ${summary} -T ${tree} -l ${split} 2>&1 | tee -a "${dest}"

        done
    done
done

current_date_time="$(date "+%d %H:%M:%S")"
echo "$current_date_time"
