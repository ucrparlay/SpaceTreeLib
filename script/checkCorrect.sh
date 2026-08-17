#!/bin/bash
#
# Compares the trees against the CGAL oracle and exits non-zero on the first
# mismatch.
#
#   ./checkCorrect.sh [fast|full] [data-root] [build-dir]
#
# fast (default) runs one small size and takes a couple of minutes; full runs
# the whole matrix and takes hours. The data root defaults to $PSI_DATA and
# then to the author's machine; generate your own with:
#
#   build/data_generator -p <root> -n 1000000 -d 2 -file_num 2 -varden 0
#   build/data_generator -p <root> -n 1000000 -d 2 -file_num 2 -varden 1
#
# The oracle needs -DCGAL=ON, and the checks it makes are asserts, so it needs
# -DDEBUG=ON to mean anything.

set -euo pipefail

tier="${1:-fast}"
data_root="${2:-${PSI_DATA:-/ssd1/zmen002/kdtree}}"
build_dir="${3:-../build}"

if [[ ${tier} == "fast" ]]; then
	nodes=(${PSI_NODES:-1000000})
	dims=(2 3)
	k=100
	files_per_dir=1
elif [[ ${tier} == "full" ]]; then
	# 10000000 and 50000000 are deliberately out: between them they are most
	# of the wall clock, and the one full sweep that ran them finished 780
	# checks with 0 failures without either being the size that caught
	# anything. Ask for them explicitly with PSI_NODES if a change warrants
	# it.
	nodes=(${PSI_NODES:-1000000 5000000 8000000})
	dims=(2 3)
	k=100
	files_per_dir=0   # every file
else
	echo "usage: $0 [fast|full] [data-root] [build-dir]" >&2
	exit 2
fi

# tree 0/1 are kd and orth, checked by kd_ccp; tree 2 is p_tree, checked by
# p_ccp. The old script assigned its tester twice, so p_ccp never ran.
declare -A tester=([0]=kd_ccp [1]=kd_ccp [2]=p_ccp)
declare -A splits=([0]="0 3" [1]="3" [2]="1 2")
declare -A tree_dims=([0]="${dims[*]}" [1]="${dims[*]}" [2]="2")

# -t is a bitmask of the update paths, -q of the query paths. Both used to be
# hardcoded to a single value, so batch updates and range queries were never
# checked at all. Bit 1 (delete) needs bit 0 (insert): on its own it would ask
# batch_delete to remove points that were never inserted.
tags=(0 1 3 4)
query_types=(1 2)

make -C "${build_dir}" kd_ccp p_ccp

pass=0
fail=0
for tree in 0 1 2; do
	for path in "${data_root}/ss_varden" "${data_root}/uniform"; do
		for node in "${nodes[@]}"; do
			for dim in ${tree_dims[${tree}]}; do
				dir="${path}/${node}_${dim}"
				if [[ ! -d ${dir} ]]; then
					echo "skip: ${dir} does not exist"
					continue
				fi
				for split in ${splits[${tree}]}; do
					for tag in "${tags[@]}"; do
						for q in "${query_types[@]}"; do
							seen=0
							for file in "${dir}"/*.in; do
								[[ -e ${file} ]] || continue
								seen=$((seen + 1))
								if [[ ${files_per_dir} -gt 0 && ${seen} -gt ${files_per_dir} ]]; then
									break
								fi
								out=$("${build_dir}/${tester[${tree}]}" \
									-p "${file}" -i 1 -s 0 \
									-d "${dim}" -k "${k}" \
									-t "${tag}" -r 2 \
									-l "${split}" -T "${tree}" \
									-q "${q}" 2>&1) || true
								if grep -qi "ok" <<<"${out}"; then
									pass=$((pass + 1))
								else
									fail=$((fail + 1))
									echo "FAIL tree=${tree} split=${split} tag=${tag} q=${q} dim=${dim} ${file}"
									echo "${out}" | tail -20
									echo "re-run: ${build_dir}/${tester[${tree}]} -p ${file} -i 1 -s 0 -d ${dim} -k ${k} -t ${tag} -r 2 -l ${split} -T ${tree} -q ${q}"
									exit 1
								fi
							done
						done
					done
				done
			done
		done
	done
done

echo "PASS=${pass} FAIL=${fail}"
if [[ ${fail} -ne 0 ]]; then
	exit 1
fi
# Reporting success having run nothing is how this script used to pass.
if [[ ${pass} -eq 0 ]]; then
	echo "no configuration ran: is ${data_root} populated?" >&2
	exit 1
fi
