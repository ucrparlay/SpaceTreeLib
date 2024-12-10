import os
import sys
import csv

print(os.getcwd())

path = "../benchmark"
benchmarks = ["real_world"]
storePrefix = "data/"
Nodes = [1000000000]
Dims = [2, 3]

solverName = ["kd", "orth" "r"]

benchmarks = {
    "osm": "2",
    "Cosmo50": "3",
    "GeiLifeNoScale": "3",
    "Household": "7",
    "HT": "10",
}

resMap = {
    "kd": "res_0_summary.out",
    "orth": "res_1_summary.out",
    "r": "res_2_summary.out",
}

common = ["solver", "file"]

#! order by test order
files = ["real_world"]

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
    "real_world": summary_header,
}

prefix = [0] * len(files)


sep_lines = []


def combine(P, file, csvWriter, solver):
    print(P)
    lines = open(P, "r").readlines()
    if len(lines) == 0:
        return
    for line in lines:
        l = " ".join(line.split())
        l = l.split(" ")
        if l[0].endswith(".in"):
            l[0] = l[0][:-3]
            continue
        sep_lines.append([solver] + l)


def write(csvWriter):
    for line in sep_lines:
        csvWriter.writerow(line)


def csvSetup(file):
    csvFilePointer = open(storePrefix + file + ".csv", "w", newline="")
    csvFilePointer.truncate()
    csvWriter = csv.writer(csvFilePointer)
    csvWriter.writerow(common + file_header[file])
    return csvWriter, csvFilePointer


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
        csvWriter, csvFilePointer = csvSetup(file)
        for bench in benchmarks:
            for solver in solverName:
                P = path + "/" + bench + "/" + resMap[solver]
                combine(P, file, csvWriter, solver)
        csvFilePointer.close()

    # reorder()
