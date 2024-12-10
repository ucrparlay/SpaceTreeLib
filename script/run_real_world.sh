#!/bin/bash
set -o xtrace

# Solvers=("test" "zdtree" "cgal")
Solvers=("test")
Tree=(0 1 2)
DataPath="/data/legacy/data3/zmen002/kdtree/geometry"
declare -A file2Dims
file2Dims["Cosmo50"]="3"
file2Dims["GeoLifeNoScale"]="3"
file2Dims["osm"]="2"
# file2Dims["Household"]="7"
# file2Dims["HT"]="10"

tag=0
k=10
onecore=0
readFile=0
summary=0
queryType=$((2#111)) # 1110000
# queryType=1 # 1110000
type="real_world"
resFile=""

for solver in ${Solvers[@]}; do
	exe="../build/${solver}"

	for tree in "${Tree[@]}"; do
		#* decide output file
		if [[ ${solver} == "rtree" ]]; then
			resFile="rtree.out"
		elif [[ ${solver} == "test" ]]; then
			resFile="res_${tree}_${type}.out"
		fi

		log_path="../benchmark/real_world"
		mkdir -p ${log_path}
		dest="${log_path}/${resFile}"
		: >${dest}
		echo ">>>${dest}"

		for filename in "${!file2Dims[@]}"; do
			echo ${filename}

			numactl -i all ${exe} -p "${DataPath}/${filename}.in" -k ${k} -t ${tag} -d ${file2Dims[${filename}]} -q ${queryType} -i ${readFile} -s ${summary} -T ${tree} | tee -a "${dest}"

		done
	done
done

current_date_time="$(date "+%d %H:%M:%S")"
echo "$current_date_time"
