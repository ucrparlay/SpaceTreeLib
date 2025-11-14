#!/bin/bash
#
set -o xtrace

DATA_PREFIX="${1:-/data/zmen002/kdtree}"
NODE_SIZE="${2:-1000000000}"
DIMENSION="${3:-2}"

make -C ../build/ data_generator data_washer

./../build/data_generator -p ${DATA_PREFIX} -n ${NODE_SIZE} -d ${DIMENSION} -file_num 2 -varden 0
./../build/data_generator -p ${DATA_PREFIX} -n ${NODE_SIZE} -d ${DIMENSION} -file_num 2 -varden 1
./../build/data_washer -coord_type 0 -d ${DIMENSION} -usage 2 -p "${DATA_PREFIX}/uniform_bigint/${NODE_SIZE}_${DIMENSION}/2.in" -output_suffix "_sort_by_0.in"

# mkdir -p "${DATA_PREFIX}/geometry"
# wget -O "${DATA_PREFIX}/geometry/Cosmo50_round_no_dup.in" "https://www.dropbox.com/scl/fi/noh6nw2xl1ymtqtqrvgsu/Cosmo50_round_no_dup.in?rlkey=vggpfsy5v2iles0agaz2fa153&st=n3v6b4c7&dl=1"
# wget -O "${DATA_PREFIX}/geometry/GeoLifeNoScale_round_no_dup.in" "https://www.dropbox.com/scl/fi/qp22hj33tjvi2iak87agp/GeoLifeNoScale_round_no_dup.in?rlkey=t10yep32bg3hvi7jyjtug2pq8&st=z39125be&dl=1"
# wget -O "${DATA_PREFIX}/geometry/osm_round_no_dup.in" "https://www.dropbox.com/scl/fi/0op87nm597d94392is9kw/osm_round_no_dup.in?rlkey=55fei6om4vs32h1pyir2digdh&st=2i147czc&dl=1"
