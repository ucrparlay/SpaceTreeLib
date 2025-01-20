#!/bin/bash
set -o xtrace

Nodes=(1000000 5000000 8000000 10000000 50000000)
# Nodes=(50000000)

K=100
tester="ccp"
resFile="Correct.out"
dest="logger.in"
out="log.in"
: >${dest}
tag=$((2#111))
# tag=$((2#0))
count=1 # count the number of ok in the output
dims=(2 3 5)
queryTypes=$((2#111))
trees=(1 0)

Paths=("/ssd1/zmen002/kdtree/ss_varden/" "/ssd1/zmen002/kdtree/uniform/")

#* check node
for queryType in "${queryTypes[@]}"; do
	for path in "${Paths[@]}"; do
		for node in "${Nodes[@]}"; do
			for dim in "${dims[@]}"; do
				for tree in "${trees[@]}"; do
					if [[ ${tree} -eq 0 ]]; then
						# splits=(0 3)
						splits=(3)
					elif [[ ${tree} -eq 1 ]]; then
						splits=(3)
					elif [[ ${tree} -eq 2 ]]; then
						splits=(0)
					fi

					for split in "${splits[@]}"; do
						if [[ ${node} -ge 50000000 ]]; then
							K=10
						fi

						files_path="${path}${node}_${dim}"
						echo "${files_path}"

						for file in "${files_path}/"*.in; do
							echo "------->${file}"
							../build/${tester} -p "${file}" -i 1 -s 0 -d "${dim}" -k ${K} -t ${tag} -r 2 -l "${split}" -T "${tree}" -q ${queryType} >>${dest}

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
done

echo "OK, Well done :)"
echo "OK, Well done :)" >>"log.in"
