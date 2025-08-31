#!/bin/bash
#
set -o xtrace

/usr/bin/drop_caches

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

origin_paths=("/data/zmen002/kdtree/ss_varden_bigint/1000000000_$((dim))/1.in" "/data/zmen002/kdtree/uniform_bigint/1000000000_${dim}/2_sort_by_0.in" "/data/zmen002/kdtree/uniform_bigint/1000000000_${dim}/2.in")
insert_paths=("/data/zmen002/kdtree/ss_varden_bigint/1000000000_$((dim))/2.in" "/data/zmen002/kdtree/uniform_bigint/1000000000_${dim}/1_sort_by_0.in" "/data/zmen002/kdtree/uniform_bigint/1000000000_${dim}/1.in")
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
