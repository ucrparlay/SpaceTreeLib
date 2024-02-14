#!/bin/bash
# {
#     sleep 210m
#     kill $$
# } &

Solvers=("test")
# Solvers=("zdtree" "test")
# Node=(100000000 1000000000)
Node=(100000000 1000000000)
# Dim=(2 3 5)
Dim=(9)
declare -A datas
datas["/data3/zmen002/kdtree/ss_varden/"]="../benchmark/ss_varden/"
datas["/data3/zmen002/kdtree/uniform/"]="../benchmark/uniform/"

tag=2
k=10
onecore=0
insNum=2
# queryType=1
queryType=$((2#1001)) # 1110000
echo $queryType
type="summary"

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
		if [ ${dim} -gt 3 ] && [ ${solver} == "zdtree" ]; then
			continue
		fi

		for dataPath in "${!datas[@]}"; do
			for node in ${Node[@]}; do
				files_path="${dataPath}${node}_${dim}"
				log_path="${datas[${dataPath}]}${node}_${dim}"
				mkdir -p ${log_path}
				dest="${log_path}/${resFile}"
				: >${dest}
				echo ">>>${dest}"

				if [[ ${dim} == 9 ]] && [[ ${solver} == "cgal" ]]; then
					rounds=1
					insNum=1
				else
					rounds=3
					insNum=2
				fi

				for ((i = 1; i <= ${insNum}; i++)); do

					export PARLAY_NUM_THREADS=192
					timeout 11200s numactl -i all ${exe} -p "${files_path}/${i}.in" -k ${k} -t ${tag} -d ${dim} -q ${queryType} -r ${rounds} >>${dest}

					retval=$?
					if [ ${retval} -eq 124 ]; then
						echo -e "timeout" >>${dest}
						echo "timeout ${node}_${dim}"
					else
						echo "finish ${node}_${dim}"
					fi
				done
			done
		done
	done
done
