#!/usr/bin/env python3
import re
import sys
import csv
import io
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass, field
from pathlib import Path

incre_ratios = [1, 0.1, 0.01, 0.001, 0.0001]
solver = ""
ratio_map = {}
input_path = sys.argv[1]
store_prefix = "data/"
type = "incre_insert" if input_path.find(
    "incre_insert") != -1 else "incre_delete"


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
        to_float(tot_num)
    ]


def parse_knn_time(text):
    """
    Extracts all numbers after each 'knn time:' label from the input text.
    Returns a dictionary {label: [list of floats]}.
    """
    if text == [""]:
        return [0.0, 0.0, 0.0] * 4

    # print(text)
    pattern = r'(\S+)\s+knn time:\s*((?:[\d\.]+\s*)+)'
    result = {}
    query_type_order = ['in-dis-skewed', 'out-dis-skewed',
                        'in-dis-uniform', 'out-dis-uniform']
    for line in text:
        match = re.search(pattern, line)
        if match:
            label = match.group(1)
            numbers_str = match.group(2)
            numbers = [float(num)
                       for num in re.findall(r'[\d\.]+', numbers_str)]
            result[label] = [numbers[i]
                             for i in range(0, len(numbers), 6)] if numbers else [0.0, 0.0, 0.0]
    sorted_res = {key: result[key]
                  for key in query_type_order if key in result}
    knn_res = list(sorted_res.values())
    flattened = [item for sublist in knn_res for item in sublist]
    return flattened


def combine(P, csv_writer):
    lines = open(P, "r").readlines()
    if len(lines) == 0:
        return

    solver = ""
    benchmark = ""
    index = 0
    while index < len(lines):
        lin_sep = " ".join(lines[index].split())
        lin_sep = lin_sep.split(" ")
        # print(lin_sep)

        if len(lin_sep) == 0 or lin_sep[0] == "":  # Skip empty lines
            index += 1
            continue

        # Tree information
        if lin_sep[0] == 'Tree:':
            if lin_sep[1] == "KdTree;":
                solver = "KdTree"
            elif lin_sep[1] == "OrthTree;":
                solver = "OrthTree"
            elif lin_sep[1] == "PTree;":
                solver = "PTree"
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
                type = lin_sep[1]
                index += 1
            case "##":
                data = parse_insert_time(lines[index + 1])
                ratio = float(lin_sep[1])
                if (lin_sep[1] == "1"):
                    data = data + parse_knn_time([""])
                    index += 2
                else:
                    data = data + parse_knn_time(lines[index+2:index+2+4])
                    index += 6
                ratio_map[ratio] = data
                csv_writer.writerow(
                    [benchmark, solver, ratio] + data
                )


def csvSetup():
    csv_file_pointer = open(store_prefix + type + ".csv", "w", newline="")
    print(csv_file_pointer.name)
    csv_file_pointer.truncate()
    csv_writer = csv.writer(csv_file_pointer)
    csv_writer.writerow(['benchmark', 'solver', 'ratio', 'median',
                         'avg', 'min', 'max', 'tot'] + ['k=1', 'k=2', 'k=3']*3)
    return csv_writer, csv_file_pointer


csv_writer, csv_file_pointer = csvSetup()
combine(input_path, csv_writer)
