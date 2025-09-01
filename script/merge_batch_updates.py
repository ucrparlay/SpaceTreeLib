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
    f"cat $(ls logs/batch_updates/* | grep -v -e {input_type}.log) > logs/batch_updates/{input_type}.log",
    shell=True,
)
input_path = "logs/batch_updates/" + input_type + ".log"
store_path = "data/batch_updates/" + input_type + ".csv"
type = ""


def to_float(match):
    return float(match.group(1)) if match else None


def combine(P) -> List:
    lines = open(P, "r").readlines()

    solver = ""
    benchmark = ""
    index = 0
    all_data = []
    op = ""
    while index < len(lines):
        lin_sep = " ".join(lines[index].split())
        lin_sep = lin_sep.split(" ")
        # print(lin_sep)

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

        # Query information
        data = [[benchmark, solver] + lines[index].split()]
        index += 1

        all_data = all_data + data
    return all_data


def post_processing(data):
    # desired order for x[0]
    bench_order = ["Uniform_by_x", "Varden", "Uniform"]
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
            x[2],
        ),
    )
    return sorted_lst


def csvSetup():
    csv_file_pointer = open(store_path, "w", newline="")
    print(csv_file_pointer.name)
    csv_file_pointer.truncate()
    csv_writer = csv.writer(csv_file_pointer)
    csv_writer.writerow(
        ["benchmark", "solver"] +
        ["i0.1", "i0.2", "i0.5", "i1", "i2", "i5", "i10",
         "i20", "i50", "i100", "i200", "i500", "i1000"]
        + ["d0.1", "d0.2", "d0.5", "d1", "d2", "d5", "d10",
           "d20", "d50", "d100", "d200", "d500", "d1000"]
    )
    return csv_writer, csv_file_pointer


csv_writer, csv_file_pointer = csvSetup()
data = combine(input_path)
# boost_res = combine_boost("logs/batch_updates/2_4_0.log")
# data = data + boost_res
data = post_processing(data)
# print(data)
csv_writer.writerows(data)
