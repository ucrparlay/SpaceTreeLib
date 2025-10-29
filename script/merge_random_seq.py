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
import glob

# Get script directory to handle relative paths
script_dir = os.path.dirname(os.path.abspath(__file__))
log_folder_name = "random_seq_pair"
log_dir = os.path.join(script_dir, "logs", log_folder_name)
output_dir = os.path.join(script_dir, "data", log_folder_name)
output_file = os.path.join(output_dir, "summary.csv")


def parse_filename(filename: str) -> Tuple[str, str, int]:
    """
    Parse filename to extract benchmark, random_order (ON/OFF) and leaf_wrap value.
    Expected format: {benchmark}_ON_{leaf_wrap}.log or {benchmark}_OFF_{leaf_wrap}.log
    Example: 1.in_ON_32.log -> ("1.in", "ON", 32)
    """
    basename = os.path.basename(filename)
    name_without_ext = basename.replace(".log", "")
    parts = name_without_ext.split("_")

    # Format: benchmark_ON/OFF_number
    # Example: 1.in_ON_32 or 2.in_OFF_1024
    if len(parts) >= 3:
        # Last part is leaf_wrap number
        leaf_wrap = int(parts[-1])
        # Second to last is ON/OFF
        random_order = parts[-2]
        # Everything before that is the benchmark
        benchmark = "_".join(parts[:-2])
        return benchmark, random_order, leaf_wrap

    return None, None, None


def extract_data_from_file(filepath: str) -> Tuple[Optional[int], Optional[float]]:
    """
    Extract leaf_wrap from first line and time from second line of the file.
    First line format: Tree: CPAM; AugType: HasBox; Split: HilbertCurve; BDO: 6; Inba: 30; LeafWrap: 32; Coord: integer;
    Second line format: 1.in 8.29854 -1 -1 24.1705
    Returns: (leaf_wrap_from_file, time)
    """
    try:
        with open(filepath, "r") as f:
            lines = f.readlines()
            if len(lines) < 2:
                print(f"Warning: File {filepath} has less than 2 lines")
                return None, None

            # Extract LeafWrap from first line
            first_line = lines[0].strip()
            leaf_wrap_match = re.search(r"LeafWrap:\s*(\d+)", first_line)
            leaf_wrap_from_file = (
                int(leaf_wrap_match.group(1)) if leaf_wrap_match else None
            )

            # Extract time from second line (second element)
            second_line = lines[1].strip()
            parts = second_line.split()
            time_value = None
            if len(parts) >= 2:
                try:
                    time_value = float(parts[-1])
                except ValueError:
                    print(
                        f"Warning: Could not parse time from: {parts[-1]} in {filepath}"
                    )

            return leaf_wrap_from_file, time_value

    except Exception as e:
        print(f"Error reading file {filepath}: {e}")
        return None, None


def process_logs() -> List[List]:
    """
    Process all log files in the log directory.
    """
    all_data = []

    # Mapping from file names to readable benchmark names
    benchmark_names = {"1.in": "Varden", "2.in": "Uniform"}

    # Mapping from ON/OFF to RAND/SORT
    order_names = {"ON": "RAND", "OFF": "SORT"}

    # Find all log files
    log_files = glob.glob(os.path.join(log_dir, "*.log"))

    if not log_files:
        print(f"Warning: No log files found in {log_dir}")
        return all_data

    print(f"Found {len(log_files)} log files in {log_dir}")

    for log_file in sorted(log_files):
        # Parse filename
        benchmark, random_order, leaf_wrap_filename = parse_filename(log_file)

        if benchmark is None:
            print(f"Warning: Could not parse filename: {log_file}")
            continue

        # Map benchmark name and order name
        benchmark_display = benchmark_names.get(benchmark, benchmark)
        order_display = order_names.get(random_order, random_order)

        # Extract data from file content
        leaf_wrap_file, time_value = extract_data_from_file(log_file)

        # Verify leaf_wrap matches
        if leaf_wrap_file is not None and leaf_wrap_filename != leaf_wrap_file:
            print(
                f"Warning: LeafWrap mismatch in {os.path.basename(log_file)}: "
                f"filename={leaf_wrap_filename}, file={leaf_wrap_file}"
            )

        print(
            f"  {os.path.basename(log_file)}: benchmark={benchmark_display}, "
            f"order={order_display}, leaf_wrap={leaf_wrap_filename}, time={time_value}"
        )

        # Add to data
        all_data.append(
            [benchmark_display, order_display, leaf_wrap_filename, time_value]
        )

    return all_data


def write_csv(data: List[List]):
    """
    Write data to CSV file.
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Sort data by benchmark, leaf_wrap (as integer), and order
    sorted_data = sorted(data, key=lambda x: (x[0], int(x[2]), x[1]))

    with open(output_file, "w", newline="") as f:
        csv_writer = csv.writer(f)

        # Write header
        csv_writer.writerow(["benchmark", "order", "leaf_wrap", "time"])

        # Write data
        csv_writer.writerows(sorted_data)

    print(f"\nCSV file written to: {output_file}")


def main():
    # Process all log files
    data = process_logs()

    if not data:
        print("No data to write")
        return

    # Write to CSV
    write_csv(data)

    print(f"Processed {len(data)} log files successfully")


if __name__ == "__main__":
    main()
