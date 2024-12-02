#!/bin/bash
set -o xtrace

# ./gen_ss_varden.sh -g 10 -n 1000000 -d 3 -v 1

# NOTE: varden
./../build/data_generator -p /data/zmen002/kdtree/ -n 1000000000 -d 2 -t 1 -f 2
./../build/data_generator -p /data/zmen002/kdtree/ -n 1000000000 -d 3 -t 1 -f 2

# NOTE: uniform
./../build/data_generator -p /data/zmen002/kdtree/ -n 1000000000 -d 2 -t 0 -f 2
./../build/data_generator -p /data/zmen002/kdtree/ -n 1000000000 -d 3 -t 0 -f 2
