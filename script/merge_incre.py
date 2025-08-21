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
input_path = "logs/incre_update/" + input_type + ".log"
store_path = "data/incre_update/" + input_type + ".csv"
type = ""


def to_float(match):
    return float(match.group(1)) if match else None


def parse_insert_time(text):
    """
    Extracts the 2nd number after median, min, max, and the number after tot and avg.
    Returns a list of [median, min, max, tot, avg] as floats.
    """

    tuple_pattern = r"\(\s*\d+\s*,\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*\)"
    median_num = re.search(r"median:\s*" + tuple_pattern, text)
    min_num = re.search(r"min:\s*" + tuple_pattern, text)
    max_num = re.search(r"max:\s*" + tuple_pattern, text)

    # Regex for single number extraction (tot/avg, covers scientific notation)
    single_pattern = r"([-\+]?\d*\.?\d+(?:[eE][-+]?\d+)?)"

    tot_num = re.search(r"tot:\s*" + single_pattern, text)
    avg_num = re.search(r"avg:\s*" + single_pattern, text)
    return [
        to_float(median_num),
        to_float(avg_num),
        to_float(min_num),
        to_float(max_num),
        to_float(tot_num),
    ]


def parse_knn_time(text):
    """
    Extracts all numbers after each 'knn time:' label from the input text.
    Returns a dictionary {label: [list of floats]}.
    """
    if text == [""]:
        return [0.0, 0.0, 0.0] * 4

    # print(text)
    pattern = r"(\S+)\s+knn time:\s*((?:[\d\.]+\s*)+)"
    result = {}
    query_type_order = [
        "in-dis-skewed",
        "out-dis-skewed",
        "in-dis-uniform",
        "out-dis-uniform",
    ]
    for line in text:
        match = re.search(pattern, line)
        if match:
            label = match.group(1)
            numbers_str = match.group(2)
            numbers = [float(num)
                       for num in re.findall(r"[\d\.]+", numbers_str)]
            sep = 6
            if len(numbers) == 3:  # for boost r-tree
                sep = 1
            result[label] = (
                [numbers[i] for i in range(0, len(numbers), sep)]
                if numbers
                else [0.0, 0.0, 0.0]
            )
    sorted_res = {key: result[key]
                  for key in query_type_order if key in result}
    knn_res = list(sorted_res.values())
    flattened = [item for sublist in knn_res for item in sublist]
    return flattened


def parse_range_time(text):
    result = []

    for line in text:
        # Skip empty lines or lines that don't start with "range"
        if not line.strip() or not line.startswith("range"):
            continue

        # Split the line into parts and extract numbers
        parts = line.split()

        # Remove the "range" prefix and get the numbers
        numbers = []
        for part in parts[2:]:  # Skip "range" and the type (count/query)
            try:
                # Try to convert to float
                numbers.append(float(part))
            except ValueError:
                continue

        # Extract 1st, 7th, and 13th elements (0-indexed: 0, 6, 12)
        if len(numbers) == 3:  # for boost r-tree
            extracted = [numbers[0], numbers[1], numbers[2]]
            result.append(extracted)

        if len(numbers) >= 13:
            extracted = [numbers[0], numbers[6], numbers[12]]
            result.append(extracted)

    flattened = [item for sublist in result for item in sublist]
    return flattened


def combine_boost(P) -> List:
    lines = open(P, "r").readlines()

    solver = ""
    benchmark = ""
    index = 0
    all_data = []
    while index < len(lines):
        lin_sep = " ".join(lines[index].split())
        lin_sep = lin_sep.split(" ")
        # print(lin_sep)

        if len(lin_sep) == 0 or lin_sep[0] == "":  # Skip empty lines
            index += 1
            continue

        # Tree information
        if lin_sep[0] == "Tree:":
            solver = "Boost"
            match lines[index + 1].split(" ")[0]:
                case "1.in":
                    benchmark = "Varden"
                case "2_sort_by_0.in":
                    benchmark = "Uniform_by_x"
                case "2.in":
                    benchmark = "Uniform"
            index += 2

        def process(data, lin_sep_arg, index_arg):
            data = data + parse_knn_time(lin_sep_arg[index_arg: index_arg + 4])
            index_arg += 4

            data = data + \
                parse_range_time(lin_sep_arg[index_arg: index_arg + 2])
            index_arg += 2
            return data

        # Query information
        match lin_sep[0]:
            case "##":
                data = []
                ratio = 0.0
                if lin_sep[3] == "full:":
                    data = [-1, -1, -1, -1, float(lin_sep[-1])]
                    ratio = float(1)
                    data = process(data, lines, index+15)
                    index += 7
                elif ("insert" in str(input_type)) and lin_sep[2] == "insert":
                    data = [-1, -1, -1, -1, float(lin_sep[-1])]
                    ratio = float(0.0001)
                    data = process(data, lines, index+1)
                    index += 7
                elif ("delete" in str(input_type)) and lin_sep[2] == "delete":
                    data = [-1, -1, -1, -1, float(lin_sep[-1])]
                    ratio = float(0.0001)
                    data = process(data, lines, index+1)
                    index += 7
                else:
                    index += 7
                    continue

                ratio_map[ratio] = data
                all_data = all_data + [[benchmark, solver, ratio] + data]

    return all_data


def combine(P) -> List:
    lines = open(P, "r").readlines()

    solver = ""
    benchmark = ""
    index = 0
    all_data = []
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
        match lin_sep[0]:
            case "#":
                # type = lin_sep[1]
                index += 1
            case "##":
                data = parse_insert_time(lines[index + 1])
                ratio = float(lin_sep[1])

                data = data + parse_knn_time(lines[index + 2: index + 2 + 4])
                index += 6

                data = data + parse_range_time(lines[index: index + 2])
                index += 2

                ratio_map[ratio] = data
                all_data = all_data + [[benchmark, solver, ratio] + data]
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
        "MVZD",
        "Boost",
    ]  # desired order for x[1]

    sorted_lst = sorted(
        data,
        key=lambda x: (
            bench_order.index(x[0]),
            tree_order.index(x[1]),
            -x[2],
        ),
    )
    return sorted_lst


def csvSetup():
    csv_file_pointer = open(store_path, "w", newline="")
    print(csv_file_pointer.name)
    csv_file_pointer.truncate()
    csv_writer = csv.writer(csv_file_pointer)
    csv_writer.writerow(
        ["benchmark", "solver", "ratio", "median", "avg", "min", "max", "tot"]
        + ["IDS1", "IDS10", "IDS100"]
        + ["ODS1", "ODS10", "ODS100"]
        + ["IDU1", "IDU10", "IDU100"]
        + ["ODU1", "ODU10", "ODU100"]
        + ["RCS", "RCM", "RCL"]
        + ["RRS", "RRM", "RRL"]
    )
    return csv_writer, csv_file_pointer


subprocess.run(
    f"cat $(ls logs/incre_update/*{input_type}* | grep -v {input_type}.log) > logs/incre_update/{input_type}.log",
    shell=True,
)

csv_writer, csv_file_pointer = csvSetup()
data = combine(input_path)
data = data + combine_boost("logs/boost_rtree/uniform_by_x.log") + combine_boost(
    "logs/boost_rtree/varden.log") + combine_boost("logs/boost_rtree/uniform.log")
# print(data)
data = post_processing(data)
# print(data)
csv_writer.writerows(data)
