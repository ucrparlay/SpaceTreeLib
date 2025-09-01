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

incre_ratios = [1, 0.1, 0.01, 0.001, 0.0001]
solver = ""
ratio_map = {}
input_type = sys.argv[1]

subprocess.run(
    f"""
    > logs/scalability/{input_type}.log
    for file in $(ls logs/scalability/*.in_*.log); do
        filename=$(basename "$file")
        number=$(echo "$filename" | grep -o '\.in_[0-9]\+' | cut -d'_' -f2)

        sed "3s/\$/ $number/" "$file" >> logs/scalability/{input_type}.log
        echo "" >> logs/scalability/{input_type}.log
    done
    """,
    shell=True,
)

input_path = "logs/scalability/" + input_type + ".log"
store_path = "data/scalability/" + input_type + ".csv"
type = ""


def to_float(match):
    return float(match.group(1)) if match else None


def combine(P) -> List:
    lines = open(P, "r").readlines()

    solver = ""
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
            if lin_sep[1] == "KdTree;":
                solver = "KdTree"
            elif lin_sep[1] == "OrthTree;":
                solver = "OrthTree"
            elif lin_sep[1] == "PTree;":
                solver = "PTree-H" if lin_sep[5] == "HilbertCurve;" else "PTree-Z"
            elif lin_sep[1] == "CPAM;":
                solver = "CPAM-H" if lin_sep[5] == "HilbertCurve;" else "CPAM-Z"
            elif lin_sep[1] == "ZdTree;":
                solver = "ZdTree"
            else:
                raise ValueError(f"Unknown solver: {lin_sep[1]}")

            match lines[index + 1].split(" ")[0]:
                case "1.in":
                    benchmark = "Varden"
                case "2_sort_by_0.in":
                    benchmark = "Uniform_by_x"
                case "2.in":
                    benchmark = "Uniform"
            index += 2
            continue
        else:
            all_data = all_data + \
                [[benchmark, solver, lin_sep[3], lin_sep[0], lin_sep[1], lin_sep[2]]]
            index += 1
    return all_data


def post_processing(data):
    # desired order for x[0]
    bench_order = ["Uniform", "Uniform_by_x", "Varden"]
    tree_order = [
        "PTree-H",
        "PTree-Z",
        "CPAM-H",
        "CPAM-Z",
        "KdTree",
        "OrthTree",
        "ZdTree",
        "Boost",
    ]  # desired order for x[1]

    sorted_lst = sorted(
        data,
        key=lambda x: (
            bench_order.index(x[0]),
            tree_order.index(x[1]),
            int(x[2])
        ),
    )
    return sorted_lst


def csvSetup():
    csv_file_pointer = open(store_path, "w", newline="")
    print(csv_file_pointer.name)
    csv_file_pointer.truncate()
    csv_writer = csv.writer(csv_file_pointer)
    csv_writer.writerow(
        ["benchmark", "solver", "threads", "build", "insert", "delete"]
    )
    return csv_writer, csv_file_pointer


csv_writer, csv_file_pointer = csvSetup()
data = combine(input_path)
print(data)
data = post_processing(data)
# print(data)
csv_writer.writerows(data)
