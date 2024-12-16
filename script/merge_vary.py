import os
import sys
import csv

print(os.getcwd())

path = "../benchmark"
# benchmarks = ["uniform", "ss_varden"]
benchmarks = ["uniform_bigint", "ss_varden_bigint"]
storePrefix = "data/"
Nodes = [1000000000]
Dims = [2, 3]

# type = "batch_update"
# type = "batch_knn_query"
type = "query"
# type = "summary"
# type = "quality"
# type = "count"

#! order by test order
files = []
solver_name = []

if type == "batch_update":
    solver_name = ["test", "zdtree", "cgal"]
    files = ["build", "insert", "delete"]
    Dims = [3]
elif type == "query":
    # solverName = ["kdtree", "orth", "r"]
    solver_name = ["kdtree", "orth"]
    files = ["build", "knn_3", "count_3", "rquery_3"]
    # files = ["build", "knn_3"]
    Dims = [2, 3]
elif type == "summary":
    solver_name = ["kdtree", "orth"]
    files = ["build", "insert", "delete", "diff", "knn_1", "count_1", "rquery_1"]
    Dims = [2, 3]
elif type == "quality":
    solver_name = ["test"]
    files = ["build", "increBuild", "decreBuild", "increKNN"]
elif type == "real_world":
    solver_name = ["test", "zdtree", "cgal", "LogTree", "BhlTree"]
    files = ["real_world"]
elif type == "count":
    solver_name = ["test"]
    files = ["build", "count"]
    Dims = [2, 3, 5, 9]

res_map = {
    "kdtree": "res_0_",
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
res_suffix = "_summary.out"

common = [
    "solver",
    "split",
    "benchType",
    "nodes",
    "dims",
]
build_header = ["build", "max-depth", "ave-depth"]
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

file_header = {
    "build": build_header,
    "insert": insert_header,
    "delete": delete_header,
    "diff": diff_header,
    "knn_1": knn_1_header,
    "count_1": count_1_header,
    "rquery_1": rquery_1_header,
    "knn_3": knn_3_header,
    "count_3": count_3_header,
    "rquery_3": rquery_3_header,
}

nodes_map = {
    100000000: "100M",
    1000000000: "1000M",
}

prefix = [0] * len(files)


def combine(P, file, csv_writer, solver, bench_name, node, dim, split):
    if not os.path.isfile(P):
        print("No file found: " + P)
        return

    lines = open(P, "r").readlines()
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
    for i in range(0, len(sep_lines)):
        if r > len(sep_lines[i]):
            print(
                f"l ({l}), r ({r}), solver ({solver}), bench ({bench_name}), node ({node}), dim ({dim})"
            )
            print(f"len(sep_lines[{i}]): {len(sep_lines[i])}")
            raise ValueError(
                f"r ({r}) is greater than the length of sep_lines[{i}] ({len(sep_lines[i])})"
            )

    # num = 2 if solverName.index(solver)<=2 else 1
    num = int(len(lines) / 2)
    for i in range(0, len(sep_lines), num):
        line = [0.0] * width
        for j in range(i, num):
            for k in range(l, r):
                line[k - l] = line[k - l] + float(sep_lines[j][k]) / num

            # print([solver, benchName, node, dim] + list(map(lambda x: round(x, 5), line)))

        csv_writer.writerow(
            [solver, split_map[split], bench_name, nodes_map[node], dim]
            + list(map(lambda x: round(x, 5), line))
        )


def csvSetup(file):
    print(file)
    csvFilePointer = open(storePrefix + file + ".csv", "w", newline="")
    csvFilePointer.truncate()
    csvWriter = csv.writer(csvFilePointer)
    csvWriter.writerow(common + file_header[file])
    return csvWriter


def calculatePrefix():
    for i in range(0, len(files), 1):
        prefix[i] = len(file_header[files[i]])
    l = prefix[0]
    prefix[0] = 1
    for i in range(1, len(files), 1):
        r = prefix[i]
        prefix[i] = l + prefix[i - 1]
        l = r


# * merge the result
if len(sys.argv) > 1 and int(sys.argv[1]) == 1:
    print(type)
    calculatePrefix()
    for file in files:
        csvWriter = csvSetup(file)

        for bench in benchmarks:
            for dim in Dims:
                for node in Nodes:
                    for solver in solver_name:
                        for split in split_map:
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
                            combine(P, file, csvWriter, solver, bench, node, dim, split)

print(">>>ok, done")
