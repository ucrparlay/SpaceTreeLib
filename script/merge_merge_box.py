#!/usr/bin/env python3
import re
import subprocess
import sys
import csv
import os
import io
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass, field
from pathlib import Path

solver = ""
ratio_map = {}
input_type = sys.argv[1]
log_folder_name = "merge_box"

subprocess.run(
    f"""
    > logs/{log_folder_name}/{input_type}.log
    for file in $(ls logs/{log_folder_name}/*.in_*.log); do
        filename=$(basename "$file")
        cat "$file" >> logs/{log_folder_name}/{input_type}.log
        echo "" >> logs/{log_folder_name}/{input_type}.log
    done
    """,
    shell=True,
)

input_path = "logs/" + log_folder_name + "/" + input_type + ".log"
store_path = "data/" + log_folder_name + "/" + input_type + ".csv"
type = ""

solver_num = 3


def parse_solvers(lines, index, solver_num, benchmark, leaf_wrap, all_data):
    for i in range(solver_num):
        lin_sep = " ".join(lines[index + i].split())
        lin_sep = lin_sep.split(" ")
        solver = lin_sep[0].strip(":")

        all_data.append([benchmark, solver, leaf_wrap] + [lin_sep[1]])


def to_float(match):
    return float(match.group(1)) if match else None


def combine(P) -> List:
    lines = open(P, "r").readlines()

    solver = "Array"
    leaf_wrap = 0
    benchmark = ""
    index = 0
    all_data = []
    while index < len(lines):
        lin_sep = " ".join(lines[index].split())
        lin_sep = lin_sep.split(" ")
        print(lin_sep)

        if len(lin_sep) == 0 or lin_sep[0] == "":  # Skip empty lines
            index += 1
            continue

        # Tree information
        if lin_sep[0] == "Tree:":
            leaf_wrap = int(lin_sep[11].strip(";"))

            match lines[index + 1].split(" ")[0]:
                case "1.in":
                    benchmark = "Varden"
                case "2_sort_by_0.in":
                    benchmark = "Uniform_by_x"
                case "2.in":
                    benchmark = "Uniform"
            index += 2

        else:
            # parse result
            parse_solvers(lines, index, solver_num, benchmark, leaf_wrap, all_data)
            index += solver_num + 1

    return all_data


def post_processing(data):
    # desired order for x[0]
    bench_order = ["Uniform", "Uniform_by_x", "Varden"]
    tree_order = ["Array-tree", "Pointer-tree", "Par-for"]  # desired order for x[1]

    sorted_lst = sorted(
        data,
        key=lambda x: (
            bench_order.index(x[0]),
            # tree_order.index(x[1]),
            int(x[2]),
        ),
    )
    return sorted_lst


def csvSetup():
    os.makedirs(os.path.dirname(store_path), exist_ok=True)
    csv_file_pointer = open(store_path, "w", newline="")
    print(csv_file_pointer.name)
    csv_file_pointer.truncate()
    csv_writer = csv.writer(csv_file_pointer)
    csv_writer.writerow(["benchmark", "solver", "leaf_wrap", "time"])
    return csv_writer, csv_file_pointer


csv_writer, csv_file_pointer = csvSetup()
data = combine(input_path)
# print(data)
data = post_processing(data)
# print(data)
csv_writer.writerows(data)
