#!/bin/bash
#
set -o xtrace

/usr/bin/drop_caches

# declare -a leaf_wrap=(512 1024)

tag=0
dim=2
k=10
summary=0
read_file=0
queryType=0

log_path="logs/merge_box"
mkdir -p "${log_path}"

# make -C ../build/ testers

origin_paths=("/data/zmen002/kdtree/ss_varden_bigint/1000000000_$((dim))/1.in" "/data/zmen002/kdtree/uniform_bigint/1000000000_${dim}/2_sort_by_0.in" "/data/zmen002/kdtree/uniform_bigint/1000000000_${dim}/2.in")
tree=0
solver="testers"
split=0

for j in "${!origin_paths[@]}"; do
    # Set leaf_wrap based on path type
    if [[ "${origin_paths[${j}]}" == *"ss_varden_bigint"* ]]; then
        declare -a leaf_wrap=(16 32 64 128 512 1024)
    else
        declare -a leaf_wrap=(2 4 16 32 64 128 512 1024)
    fi

    for i in "${!leaf_wrap[@]}"; do
        cur_path=$(pwd)
        cd ../build
        cmake -DDEBUG=OFF -DCGAL=ON -DCMAKE_CXX_COMPILER=g++-14 -DLEAF="${leaf_wrap[${i}]}" "~/SpaceTreeLib"
        make testers
        cd ${cur_path}

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
