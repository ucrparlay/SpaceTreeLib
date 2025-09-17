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

subprocess.run(
    f"""
    > logs/linear_wrap/{input_type}.log
    for file in $(ls logs/linear_wrap/*.in_*.log); do
        filename=$(basename "$file")
        cat "$file" >> logs/linear_wrap/{input_type}.log
        echo "" >> logs/linear_wrap/{input_type}.log
    done
    """,
    shell=True,
)

input_path = "logs/linear_wrap/" + input_type + ".log"
store_path = "data/linear_wrap/" + input_type + ".csv"
type = ""


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
            leaf_wrap = int(lin_sep[11].strip(';'))

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
            solver = "Pointer" if solver == "Array" else "Array"
            all_data.append([benchmark, solver, leaf_wrap] +
                            [lin_sep[i] for i in range(0, len(lin_sep), 6)])
            # all_data = all_data + \
            #     [[benchmark, solver, leaf_wrap, lin_sep[0],
            #         lin_sep[6], lin_sep[12], lin_sep[18], lin_sep[24], lin_sep[30], lin_sep[36]]]
            index += 1

    return all_data


def post_processing(data):
    # desired order for x[0]
    bench_order = ["Uniform", "Uniform_by_x", "Varden"]
    tree_order = [
        "Pointer",
        "Array"
    ]  # desired order for x[1]

    sorted_lst = sorted(
        data,
        key=lambda x: (
            bench_order.index(x[0]),
            # tree_order.index(x[1]),
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
        ["benchmark", "solver", "leaf_wrap", "InD1",
            "InD10", "InD100", "InD1000", "OOD1", "OOD10", "OOD100", "OOD1000", "S", "M", "L"]
    )
    return csv_writer, csv_file_pointer


csv_writer, csv_file_pointer = csvSetup()
data = combine(input_path)
# print(data)
# data = post_processing(data)
# print(data)
csv_writer.writerows(data)
