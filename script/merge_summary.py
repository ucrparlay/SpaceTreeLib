import os
import sys
import csv

print(os.getcwd())

path = "../benchmark"
# benchmarks = ["uniform", "ss_varden"]
benchmarks = ["uniform_bigint", "ss_varden_bigint"]
store_prefix = "data/"
Nodes = [1000000000]
Dims = [2, 3]

solver_name = ["kd", "orth"]
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
type = "summary"

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
                    for split in list(split_name):
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
                                + type
                                + "_"
                                + split
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
