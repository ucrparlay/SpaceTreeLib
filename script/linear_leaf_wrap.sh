#!/bin/bash
#
set -o xtrace

/usr/bin/drop_caches

declare -a leaf_wrap=(4 16 32 64 128 512 1024)

tag=0
dim=2
k=10
summary=0
read_file=0
queryType=0

log_path="logs/linear_wrap"
mkdir -p "${log_path}"

# make -C ../build/ testers

origin_paths=("/data/zmen002/kdtree/ss_varden_bigint/1000000000_$((dim))/1.in" "/data/zmen002/kdtree/uniform_bigint/1000000000_${dim}/2_sort_by_0.in" "/data/zmen002/kdtree/uniform_bigint/1000000000_${dim}/2.in")
tree=0
solver="testers"
split=0

for i in "${!leaf_wrap[@]}"; do
    cur_path=$(pwd)
    cd ../build
    cmake -DDEBUG=OFF -DCGAL=ON -DCMAKE_CXX_COMPILER=g++-14 -DLEAF="${leaf_wrap[${i}]}" "~/SpaceTreeLib"
    make testers
    cd ${cur_path}

    for j in "${!origin_paths[@]}"; do
        file="${origin_paths[${j}]##*/}"
        dest="${log_path}/${file}_${leaf_wrap[${i}]}.log"
        : >"${dest}"
        echo ">>>${dest}"
        exe="../build/${solver}"

        input="${origin_paths[${j}]}"
        echo ${input}

        ${commands[${i}]} "${exe}" -p "${input}" -k ${k} -t ${tag} -d 2 -r 3 -q 0 -i 0 -s 0 -T ${tree} -l ${split} 2>&1 | tee -a "${dest}"
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
