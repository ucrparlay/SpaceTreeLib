#!/bin/bash

Solvers=("test")
Node=(100000000)
Dim=(3)
declare -A datas
datas["/data3/zmen002/kdtree/ss_varden/"]="../benchmark/ss_varden/"
datas["/data3/zmen002/kdtree/uniform/"]="../benchmark/uniform/"

tag=0
k=100
onecore=0
insNum=2
queryType=$((2#110000)) # 1110000
type="batch_update"
resFile=""

for solver in ${Solvers[@]}; do
	exe="../build/${solver}"

	#* decide output file
	if [[ ${solver} == "test" ]]; then
		resFile="res_${type}.out"
	elif [[ ${solver} == "cgal" ]]; then
		resFile="cgal_${type}.out"
	elif [[ ${solver} == "zdtree" ]]; then
		resFile="zdtree_${type}.out"
		exe="/home/zmen002/pbbsbench_x/build/zdtree"
	fi

	for dim in ${Dim[@]}; do
		for dataPath in "${!datas[@]}"; do
			for node in ${Node[@]}; do
				files_path="${dataPath}${node}_${dim}"
				log_path="${datas[${dataPath}]}${node}_${dim}"
				mkdir -p ${log_path}
				dest="${log_path}/${resFile}"
				: >${dest}
				echo ">>>${dest}"

				for ((i = 1; i <= ${insNum}; i++)); do
					PARLAY_NUM_THREADS=192 numactl -i all ${exe} -p "${files_path}/${i}.in" -k ${k} -t ${tag} -d ${dim} -q ${queryType} -i 1 >>${dest}
				done
			done
		done
	done
done

echo date +%d%H%M%S
