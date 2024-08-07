#!/bin/bash

Nodes=(1000000 5000000 8000000 10000000 50000000)
# Nodes=(10000000 50000000)

K=100
tester="checkCorrectParallel"
resFile="Correct.out"
dest="logger.in"
out="log.in"
: >${dest}
tag=0
count=1 # count the number of ok in the output
dims=(2)
queryTypes=(0)
trees=(0 1)
# queryTypes=(0 1 2)

# Paths=("/ssd0/zmen002/kdtree/uniform_bigint/" "/ssd0/zmen002/kdtree/ss_varden/")
Paths=("/localdata/zmen002/kdtree/ss_varden/" "/localdata/zmen002/kdtree/uniform_bigint/")
# Paths=("/localdata/0/zmen002/kdtree/ss_varden/" "/localdata/0/zmen002/kdtree/uniform_bigint/")

#* check node
for queryType in "${queryTypes[@]}"; do
	for path in "${Paths[@]}"; do
		for node in "${Nodes[@]}"; do
			for dim in "${dims[@]}"; do
				for tree in "${trees[@]}"; do
					if [ "${queryType}" -gt 0 ] && [ "${node}" -gt 8000000 ]; then
						continue
					fi

					if [[ ${node} -ge 50000000 ]]; then
						K=10
					fi

					files_path="${path}${node}_${dim}"
					echo "${files_path}"

					for file in "${files_path}/"*.in; do
						echo "------->${file}"
						../build/${tester} -p ${file} -d ${dim} -k ${K} -t ${tag} -r 1 -T ${tree} -q ${queryType} >>${dest}

						nc=$(grep -i -o "ok" ${dest} | wc -l)
						if [[ ${nc} -ne ${count} ]]; then
							echo 'wrong'
							exit
						fi
						count=$((count + 1))
					done
				done
			done
		done
	done
done

echo "OK, Well done :)"
echo "OK, Well done :)" >>"log.in"
