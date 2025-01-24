import os
import sys
import csv

print(os.getcwd())

path = "../benchmark"
storePrefix = "data/real_world/"
Nodes = [1000000000]
Dims = [2, 3]

solverName = ["kd", "orth"]
res_map = {
    "kd": "res_0_",
    "orth": "res_1_",
    "r": "res_2_",
}
split_name = {"0", "1", "2", "3"}
split_map = {
    "0": "MaxStr/Obj",
    "1": "Rot/Obj",
    "2": "MaxStr/SpaMid",
    "3": "Rot/SpaMid",
}
bench_node = {
    "osm": "1298M",
    "Cosmo50": "321M",
    "GeoLifeNoScale": "24M",
    "Household": "2.04M",
    "HT": "928K",
}
bench_dim = {
    "osm": "2",
    "Cosmo50": "3",
    "GeoLifeNoScale": "3",
    "Household": "7",
    "HT": "10",
}
type = "real_world"

common = ["bench", "nodes", "dims", "solver", "split"]

#! order by test order
files = ["build", "knn_3", "count_3", "rquery_3"]

build_header = [
    "build",
    "maxDepth",
    "aveDepth",
]
knn_1_header = [
    "k=10",
    "vis",
    "gen",
    "check",
    "skip",
]
count_1_header = [
    "L",
    "vis",
    "gen",
    "full",
    "skip",
]
rquery_1_header = count_1_header
knn_3_header = [
    metric for _ in range(3) for metric in ["time", "vis", "gen", "check", "skip"]
]
count_3_header = [
    metric for _ in range(3) for metric in ["time", "vis", "gen", "full", "skip"]
]
rquery_3_header = count_3_header

summary_header = build_header + knn_3_header + count_3_header + rquery_3_header
file_header = {
    "build": build_header,
    "knn_3": knn_3_header,
    "count_3": count_3_header,
    "rquery_3": rquery_3_header,
}

prefix = [0] * len(files)


sep_lines = []


def combine(P, file, csv_writer, solver, split):
    if not os.path.isfile(P):
        print("No file fonund: " + P)
        return

    print(P)
    lines = open(P, "r").readlines()
    if len(lines) == 0:
        return
    for index, line in enumerate(lines):
        if index % 2 == 0:
            continue
        l = " ".join(line.split())
        l = l.split(" ")
        if l[0].endswith(".in"):
            l[0] = l[0][:-3]
        width = len(file_header[file])
        left = prefix[files.index(file)]
        right = left + width
        sep_lines.append(
            [l[0], bench_node[l[0]], bench_dim[l[0]]]
            + [solver]
            + [split_map[split]]
            + l[left:right]
        )


def write(csv_writer):
    sep_lines.sort(key=lambda x: x[0])
    for line in sep_lines:
        csv_writer.writerow(line)


def csvSetup(file):
    csvFilePointer = open(storePrefix + file + ".csv", "w", newline="")
    csvFilePointer.truncate()
    csv_writer = csv.writer(csvFilePointer)
    csv_writer.writerow(common + file_header[file])
    return csv_writer, csvFilePointer


def calculatePrefix():
    for i in range(0, len(files), 1):
        prefix[i] = len(file_header[files[i]])
    l = prefix[0]
    prefix[0] = 1
    for i in range(1, len(files), 1):
        r = prefix[i]
        prefix[i] = l + prefix[i - 1]
        l = r

    print(prefix)


# * merge the result
if len(sys.argv) > 1 and int(sys.argv[1]) == 1:
    calculatePrefix()
    for file in files:
        sep_lines = []
        csv_writer, csvFilePointer = csvSetup(file)
        for solver in solverName:
            for split in list(split_name):
                P = (
                    path
                    + "/real_world/"
                    + res_map[solver]
                    + type
                    + "_"
                    + split
                    + ".out"
                )
                combine(P, file, csv_writer, solver, split)
        write(csv_writer)
        csvFilePointer.close()

    # reorder()
