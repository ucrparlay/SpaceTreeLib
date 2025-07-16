import os
import sys
import csv

print(os.getcwd())

path = "../benchmark"
# benchmarks = ["uniform", "ss_varden"]
benchmarks = ["uniform_bigint", "ss_varden_bigint"]
store_prefix = "data/"
Nodes = [1000000000]
# Dims = [2, 3]
Dims = [2]

solver_name = ["kdtree", "orth", "cpam"]
res_map = {
    "kdtree": "res_0_",
    "orth": "res_1_",
    "cpam": "res_2_",
}

solver_split_map = {
    "kdtree": ["kd_0", "kd_3"],
    "orth": ["orth_3"],
    "cpam": ["cpam_1", "cpam_2"],
}

split_map = {
    "kd_0": "MaxStr/Obj",
    "kd_1": "Rot/Obj",
    "kd_2": "MaxStr/SpaMid",
    "kd_3": "Rot/SpaMid",
    "orth_0": "MaxStr/Obj",
    "orth_1": "Rot/Obj",
    "orth_2": "MaxStr/SpaMid",
    "orth_3": "Rot/SpaMid",
    "cpam_1": "Hilbert",
    "cpam_2": "Morton",
}
op = "summary"

common = [
    "solver",
    "split",
    "benchType",
    "nodes",
    "dims",
]

#! order by test order
files = ["summary"]

build_header = [
    "build",
    "maxDepth",
    "aveDepth",
]
insert_header = [
    "0.01%",
    "0.1%",
    "1%",
    "10%",
]
delete_header = [
    "0.01%",
    "0.1%",
    "1%",
    "10%",
]
diff_header = [
    "10%/1%",
    "20%/1%",
    "50%/1%",
    "100%/1%",
]
knn_header = ["k=10", "vis", "gen", "check", "skip"]
range_count = [
    "L",
    "vis",
    "gen",
    "full",
    "skip",
]
range_query = range_count
summary_header = (
    build_header
    + insert_header
    + delete_header
    + diff_header
    + knn_header
    + range_query
)
file_header = {
    "summary": summary_header,
}
nodes_map = {
    100000000: "100M",
    1000000000: "1000M",
}

prefix = [0] * len(files)


def combine(P, file, csv_writer, solver, bench_name, node, dim, split):
    if not os.path.isfile(P):
        print("No file fonund: " + P)
        if solver == "zdtree":
            csv_writer.writerow(
                [solver, bench_name, node, dim] + ["-"] * len(build_header)
            )
        elif solver == "cgal":
            csv_writer.writerow(
                [solver, bench_name, node, dim] + ["T"] * len(build_header)
            )
        return

    # NOTE: begin to read the file
    print(P)
    lines = open(P, "r").readlines()
    if len(lines) == 0:
        return
    sep_lines = []
    for index, line in enumerate(lines):
        if index % 2 == 0:
            continue
        l = " ".join(line.split())
        l = l.split(" ")
        sep_lines.append(l)

    width = len(file_header[file])
    l = prefix[files.index(file)]
    r = l + width
    num = int(len(lines) / 2)
    for i in range(0, len(sep_lines), num):
        line = [0.0] * width
        # print(sep_lines[i])
        for j in range(i, num):
            for k in range(l, r):
                # print(f"j: {j}, k: {k}, l: {l}, r: {r}")
                line[k - l] = line[k - l] + float(sep_lines[j][k]) / num

        csv_writer.writerow(
            [solver, split_map[split], bench_name, node, dim]
            + list(map(lambda x: round(x, 5), line))
        )


def csvSetup(file):
    csvFilePointer = open(store_prefix + file + ".csv", "w", newline="")
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
            for dim in Dims:
                for solver in solver_name:
                    for split in solver_split_map[solver]:
                        for node in Nodes:
                            P = (
                                path
                                + "/"
                                + bench
                                + "/"
                                + str(node)
                                + "_"
                                + str(dim)
                                + "/"
                                + res_map[solver]
                                + op
                                + "_"
                                + split[-1]
                                + ".out"
                            )
                            combine(
                                P,
                                file,
                                csvWriter,
                                solver,
                                bench,
                                nodes_map[node],
                                dim,
                                split,
                            )
        csvFilePointer.close()

print(">>>ok, done")
# reorder()
