#!/usr/bin/env python

import os
import sys
import numpy as np

# utils
REV = {"A":"T","T":"A","G":"C","C":"G","N":"N"}
def revcomp(x):
    return [REV[c] for c in x[::-1]]
def frac(x):
    return x["bad"] / (x["good"] + x["bad"])
def read_condition_dict(path):
    condition_dict = dict()
    index_list = list()
    condition_names = list()
    with open(path, "rt") as f:
        next(f)
        for line in f:
            condition, forward, reverse = line.strip().split(",")
            condition_dict[condition] = (forward, reverse)
            index_list.append(forward + reverse)
            condition_names.append(condition)
    return condition_dict, index_list, condition_names
def read_sequence(path):
    result = ""
    arrowcount = 0
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                arrowcount += 1
                if arrowcount > 1:
                    raise ValueError("File contains more than one sequence.")
                continue
            result += line
    return result
def run_pear(data_path, name, fw, rv):
    os.system(
        f"pear -f {data_path}/{fw} -r {data_path}/{rv} -o {data_path}/{name} -j 16")
def run_all_pear(data_path, read_path):
    with open(read_path) as f:
        for line in f:
            name, fw, rv = line.strip().split(",")
            run_pear(data_path, name, fw, rv)

def samkey(condition, condition_dict):
    forward, reverse = condition_dict[condition]
    return forward + reverse
def revkey(key, index_list, condition_names):
    index = index_list.index(key)
    return condition_names[index]

data_path = sys.argv[1]
sequence = read_sequence(f"{data_path}/sequence.fa")
read_path = f"{data_path}/readfiles.csv"
condition_path = f"{data_path}/barcodes.csv"
out_path = f"{data_path}/outputs.csv"
count_mut = False

condition_dict, index_list, condition_names = read_condition_dict(condition_path)
example_fw, example_rv = condition_dict[condition_names[0]]
barcode_length = len(example_fw)
run_all_pear(data_path, read_path)

replicates = []
replicate_names = []
sequence_np = np.array([c for c in sequence])
for replicate_name in os.listdir(data_path):
    if ".assembled." not in replicate_name:
        continue
    replicate = f"{data_path}/{replicate_name}"
    replicate_basename = replicate_name.strip().split("/")[0].split(".")[0]
    replicate_names.append(replicate_basename)
    bad_index_count = 0
    total_count = 0
    indexed_counts = {}
    with open(replicate) as f:
        for idx, line in enumerate(f):
            total_count += 1
            if idx % 4 == 1:
                seq = line.strip()
                index_start = seq[:barcode_length]
                index_end = seq[-barcode_length:]
                index = index_start + index_end
                if index not in index_list:
                    if revcomp(index) not in index_list:
                        bad_index_count += 1
                        continue
                seq = seq[barcode_length:-barcode_length]
                mismatch = (len(seq) - len(sequence))
                if index not in indexed_counts:
                    indexed_counts[index] = dict(good=0, bad=0)
                if mismatch != 0:
                    indexed_counts[index]["bad"] += 1
                else:
                    seq_np = np.array([c for c in seq])
                    mut = (seq_np != sequence_np).any()
                    if count_mut and mut:
                        indexed_counts[index]["bad"] += 1
                    else:
                        indexed_counts[index]["good"] += 1
    replicates.append(indexed_counts)
fracrep = [[frac(r[key]) for i, r in enumerate(replicates)] for key in index_list]
fracmean = [np.array([frac(r[key]) for r in replicates]).mean() for key in index_list]
fracstd = [np.array([frac(r[key]) for r in replicates]).std() for key in index_list]

with open(out_path, "wt") as f:
    f.write(
        "condition name,index," + ",".join(replicate_names) + "\n")
    for index, condition, reps in zip(index_list, condition_names, fracrep):
        f.write(f"{condition},{index}," + ",".join([str(r) for r in reps]) + "\n")
