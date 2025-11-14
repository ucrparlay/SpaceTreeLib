#!/bin/bash
#
set -o xtrace

DATA_PREFIX="${1:-/data/zmen002/kdtree}"
NODE_SIZE="${2:-1000000000}"
DIMENSION="${3:-2}"

# /usr/bin/drop_caches

declare -a cores=(1 4 16 32 64 112 224)
declare -a threads=(
    "PARLAY_NUM_THREADS=1"
    "PARLAY_NUM_THREADS=4"
    "PARLAY_NUM_THREADS=16"
    "PARLAY_NUM_THREADS=32"
    "PARLAY_NUM_THREADS=64"
    "PARLAY_NUM_THREADS=112"
    "PARLAY_NUM_THREADS=224"
)

declare -a commands=(
    "taskset -c 0-3:4"                 # 1 thread
    "taskset -c 0-15:4 numactl -i all" # 4 threads
    "taskset -c 0-63:4 numactl -i all" # 16 threads
    "taskset -c 0-63:2 numactl -i all" # 32 threads
    "taskset -c 0-63   numactl -i all" # 64 threads
    "taskset -c 0-111  numactl -i all" # 112 threads
    "numactl -i all"                   # 224 threads
)

tag=4
dim="${DIMENSION}"
k=10
summary=0
read_file=0
queryType=0

log_path="logs/scalability"
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
        splits=(1)
    elif [[ ${tree} -eq 3 ]]; then
        solver="baselines"
        splits=(1)
    elif [[ ${tree} -eq 5 ]]; then
        solver="baselines"
        splits=(2)
    else
        echo "ERROR: Invalid tree type ${tree}"
        exit 1
    fi

    for split in "${splits[@]}"; do
        for i in "${!threads[@]}"; do
            for j in "${!origin_paths[@]}"; do
                file="${origin_paths[${j}]##*/}"
                dest="${log_path}/${file}_${cores[${i}]}_${tree}_${split}.log"
                : >"${dest}"
                echo ">>>${dest}"
                exe="../build/${solver}"
                export "${threads[${i}]}"
                if [[ ${i} -eq 0 || ${i} -eq 1 ]]; then
                    rounds=1
                else
                    rounds=2
                fi

                input="${origin_paths[${j}]}"
                insert="${insert_paths[${j}]}"
                echo ${input}
                echo ${insert}

                ${commands[${i}]} "${exe}" -p "${input}" -I "${insert}" -k ${k} -t ${tag} -d 2 -r ${rounds} -q 0 -i 1 -s 0 -T ${tree} -l ${split} 2>&1 | tee -a "${dest}"
            done
        done
    done
done

#for solver in "${Solvers[@]}"; do
#    exe="../build/${solver}"
#    #* decide output file
#    if [[ ${solver} == "test" ]]; then
#        resFile="res.out"
#    elif [[ ${solver} == "cgal" ]]; then
#        resFile="cgal.out"
#    elif [[ ${solver} == "zdtree" ]]; then
#        resFile="zdtree.out"
#        exe="/home/zmen002/pbbsbench_x/build/zdtree"
#    fi

#    for i in "${!threads[@]}"; do
#        for dataPath in "${!datas[@]}"; do
#            for node in "${Node[@]}"; do
#                logPath=${datas[${dataPath}]}

#                files_path="${dataPath}${node}_${dim}"
#                log_path="${logPath}${node}_${dim}"
#                mkdir -p "${log_path}"
#                dest="${log_path}/${cores[${i}]}_${resFile}"
#                : >"${dest}"
#                echo ">>>${dest}"

#                export "${threads[${i}]}"
#                ${commands[${i}]} "${exe}" -p "${files_path}/1.in" -k ${k} -t ${tag} -d ${dim} -r 1 -q 0 -i 1 >>"${dest}"

#            done
#        done
#    done
#done
