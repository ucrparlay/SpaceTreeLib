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
    f"cat $(ls logs/range_query_log/* | grep -v -e {input_type}.log -e 2_4_0.log) > logs/range_query_log/{input_type}.log",
    shell=True,
)
input_path = "logs/range_query_log/" + input_type + ".log"
store_path = "data/range_query_log/" + input_type + ".csv"
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

        # Query information
        match lin_sep[0]:
            case "##":
                index += 1
                data = [
                    [benchmark, solver] + line.split()
                    for line in lines[index: index + 300]
                ]
                index += 300

                all_data = all_data + data

    return all_data


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
        match lin_sep[0]:
            case "#":
                op = lin_sep[1]
                index += 1
            case "##":
                _ = parse_insert_time(lines[index + 1])
                _ = float(lin_sep[1])
                index += 3

                # Prepend benchmark and solver to each line and store in data
                data = [
                    [benchmark, solver] + line.split()
                    for line in lines[index: index + 300]
                ]
                index += 300

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
        ["benchmark", "solver", "type", "size", "time"]
    )
    return csv_writer, csv_file_pointer


csv_writer, csv_file_pointer = csvSetup()
data = combine(input_path)
boost_res = combine_boost("logs/range_query_log/2_4_0.log")
data = data + boost_res
data = post_processing(data)
# print(data)
csv_writer.writerows(data)
