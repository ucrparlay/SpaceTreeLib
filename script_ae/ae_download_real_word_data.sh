#!/bin/bash
#
set -o xtrace

DATA_PREFIX="${1:-/data/zmen002/kdtree}"

mkdir -p "${DATA_PREFIX}/geometry"
wget -O "${DATA_PREFIX}/geometry/Cosmo50_round_no_dup.in" "https://www.dropbox.com/scl/fi/noh6nw2xl1ymtqtqrvgsu/Cosmo50_round_no_dup.in?rlkey=vggpfsy5v2iles0agaz2fa153&st=n3v6b4c7&dl=1"
wget -O "${DATA_PREFIX}/geometry/osm_round_no_dup.in" "https://www.dropbox.com/scl/fi/0op87nm597d94392is9kw/osm_round_no_dup.in?rlkey=55fei6om4vs32h1pyir2digdh&st=2i147czc&dl=1"
