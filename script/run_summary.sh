#!/bin/bash
set -o xtrace

Solvers=("test")
# Solvers=("rtree")
Node=(1000000000)
# Tree=(1)
Tree=(0 1)
Dim=(2 3)
declare -A datas
# datas["/data/legacy/data3/zmen002/kdtree/ss_varden/"]="../benchmark/ss_varden/"
# datas["/data/legacy/data3/zmen002/kdtree/uniform/"]="../benchmark/uniform/"
datas["/data/zmen002/kdtree/ss_varden_bigint/"]="../benchmark/ss_varden_bigint/"
datas["/data/zmen002/kdtree/uniform_bigint/"]="../benchmark/uniform_bigint/"

tag=$((2#111)) # 1110000
k=10
insNum=2
summary=1
read_file=1
queryType=$((2#101)) # 1110000
type="summary"
round=3
resFile=""

for solver in "${Solvers[@]}"; do
    exe="../build/${solver}"

    for tree in "${Tree[@]}"; do
        if [[ ${tree} -eq 0 ]]; then
            splits=(0 3)
        elif [[ ${tree} -eq 1 ]]; then
            splits=(3)
        elif [[ ${tree} -eq 2 ]]; then
            splits=(0)
        fi

        for split in "${splits[@]}"; do

            if [[ ${solver} == "rtree" ]]; then
                resFile="rtree.out"
            elif [[ ${solver} == "test" ]]; then
                resFile="res_${tree}_${type}_${split}.out"
            fi
            for dim in "${Dim[@]}"; do
                for dataPath in "${!datas[@]}"; do
                    for node in "${Node[@]}"; do
                        files_path="${dataPath}${node}_${dim}"
                        log_path="${datas[${dataPath}]}${node}_${dim}"
                        mkdir -p "${log_path}"
                        dest="${log_path}/${resFile}"
                        : >"${dest}"
                        echo ">>>${dest}"

                        for ((i = 1; i <= insNum; i++)); do
                            numactl -i all ${exe} -p "${files_path}/${i}.in" -r ${round} -k ${k} -i ${read_file} -s ${summary} -t ${tag} -d ${dim} -q ${queryType} -T ${tree} -l ${split} 2>&1 | tee -a "${dest}"
                        done
                    done
                done
            done
        done
    done
done

current_date_time="$(date "+%d %H:%M:%S")"
echo $current_date_time
