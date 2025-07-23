#!/bin/bash
set -o xtrace

/usr/bin/drop_caches

Node=(1000000000)
# Tree=(2)
Tree=(0 1 2)
Dim=(2)
declare -A datas
# datas["/data/legacy/data3/zmen002/kdtree/ss_varden/"]="../benchmark/ss_varden/"
# datas["/data/legacy/data3/zmen002/kdtree/uniform/"]="../benchmark/uniform/"
datas["/data/zmen002/kdtree/ss_varden_bigint/"]="../benchmark/ss_varden_bigint/"
datas["/data/zmen002/kdtree/uniform_bigint/"]="../benchmark/uniform_bigint/"

tag=$((2#000)) # 1110000
k=10
insNum=2
summary=0
read_file=0
queryType=$((2#111)) # 1110000
# queryType=$((2#1)) # 1110000
type="query"
round=4
resFile=""

make -C ../build -j
for tree in "${Tree[@]}"; do
    if [[ ${tree} -eq 0 ]]; then
        solver="kd_test"
        splits=(0 3)
    elif [[ ${tree} -eq 1 ]]; then
        solver="kd_test"
        splits=(3)
    elif [[ ${tree} -eq 2 ]]; then
        solver="p_test"
        splits=(1 2)
    fi
    exe="../build/${solver}"

    for split in "${splits[@]}"; do
        resFile="res_${tree}_${type}_${split}.out"

        for dim in "${Dim[@]}"; do
            for dataPath in "${!datas[@]}"; do
                for node in "${Node[@]}"; do

                    # Hilbert curve is not supported for 3D with p_test
                    # if [ ${dim} -eq 3 ] && [ ${tree} -eq 2 ] && [ ${split} -eq 1 ]; then
                    #     continue
                    # fi

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

current_date_time="$(date "+%d %H:%M:%S")"
echo $current_date_time
