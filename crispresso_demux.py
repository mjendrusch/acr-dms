#!/usr/bin/env python

import os
import sys
import numpy as np

# utils
REV = {"A":"T","T":"A","G":"C","C":"G","N":"N"}
def revcomp(x):
    return "".join([REV[c] for c in x[::-1]])
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
def samkey(condition, condition_dict):
    forward, reverse = condition_dict[condition]
    return forward + reverse
def revkey(key, index_list, condition_names):
    index = index_list.index(key)
    return condition_names[index]
import gzip
def run_crispresso_demux(
        condition_dict, index_list: list, condition_names: list,
        data_path, name, fw, rv,
        window=None):
    condition_items = {cond: [] for cond in condition_dict}
    bad_index = set()
    # demux fastq
    print("Demultiplexing FASTQ files...")
    with gzip.open(f"{data_path}/{fw}", "rt") as ffw,\
         gzip.open(f"{data_path}/{rv}", "rt") as frv:
        lost_index = 0
        total_reads = 0
        while True:
            try:
                # read in information
                item = dict()
                item["header_fw"] = next(ffw)
                item["sequence_fw"] = next(ffw)
                item["sep_fw"] = next(ffw)
                item["quality_fw"] = next(ffw)
                item["header_rv"] = next(frv)
                item["sequence_rv"] = next(frv)
                item["sep_rv"] = next(frv)
                item["quality_rv"] = next(frv)
                total_reads += 1
                barcode_fw = item["sequence_fw"].strip()[:barcode_length]
                barcode_rv = revcomp(item["sequence_rv"].strip()[:barcode_length])
                index = barcode_fw + barcode_rv
                index_base = index
                index_rev = revcomp(index)
                if index not in index_list:
                    index = index_rev
                if index not in index_list:
                    lost_index += 1
                    bad_index.add(index_base)
                    bad_index.add(index_rev)
                    continue
                condition = condition_names[index_list.index(index)]
                condition_items[condition].append(item)
            except StopIteration:
                break
    print("Demultiplexed FASTQ files.")
    print(f"# total reads: {total_reads}")
    print(f"# lost read pairs due to invalid index: {lost_index}")
    print(f"{100 * lost_index / total_reads:.1f} % of reads lost.")
    # process fastq files
    condition_results = dict()
    os.makedirs(f"{data_path}/demultiplexed_fastq/", exist_ok=True)
    os.makedirs(f"{data_path}/CRISPResso_results/", exist_ok=True)
    for condition, items in condition_items.items():
        print(f"Running CRISPResso on {condition}")
        # write separate fastq files per condition
        r1_path = f"{data_path}/demultiplexed_fastq/{name}_{condition}_R1.fq.gz"
        r2_path = f"{data_path}/demultiplexed_fastq/{name}_{condition}_R2.fq.gz"
        crispresso_out_path = f"{data_path}/CRISPResso_results/CRISPResso_{name}_{condition}"
        with gzip.open(r1_path, "wt") as ffw,\
             gzip.open(r2_path, "wt") as frv:
            for item in items:
                ffw.write(item["header_fw"])
                ffw.write(item["sequence_fw"])
                ffw.write(item["sep_fw"])
                ffw.write(item["quality_fw"])
                frv.write(item["header_rv"])
                frv.write(item["sequence_rv"])
                frv.write(item["sep_rv"])
                frv.write(item["quality_rv"])
        # run crispresso:
        run_crispresso(r1_path, r2_path, crispresso_out_path, window=window)
        # digest crispresso:
        condition_results[condition] = crispresso_digest(
            crispresso_out_path, name, condition)
    return condition_results

def run_crispresso(r1_path, r2_path, crispresso_out_path, window=None):
    if window is None:
        window = 1
    os.system(
        f"CRISPResso -r1 {r1_path} -r2 "
        f"{r2_path} -a {sequence} -g {guide} "
        f"--cleavage_offset 0 -o {crispresso_out_path} "
        # f"--ignore_substitutions True "
        f"-p 16 "
        f"-w {window} > /dev/null 2>&1")

def crispresso_digest(crispresso_out_path, name, condition):
    with open(f"{crispresso_out_path}/CRISPResso_on_{name}_{condition}_R1_{name}_{condition}_R2/"
              f"Alleles_frequency_table_around_sgRNA_{guide}.txt", "rt") as f:
        next(f)
        p_edited = 0
        for line in f:
            algn, ref, unedited, n_del, n_ins, n_mut, n_reads, p_reads = line.strip().split("\t")
            edited = unedited != "True"
            p_reads = float(p_reads)
            if edited:
                p_edited += p_reads / 100
    return p_edited

def run_all_demux(condition_dict, index_list, condition_names, data_path, read_path,
                  window=None):
    data_by_name = dict()
    run_names = []
    with open(read_path) as f:
        for line in f:
            name, fw, rv = line.strip().split(",")
            run_names.append(name)
            p_edited = run_crispresso_demux(
                condition_dict, index_list, condition_names,
                data_path, name, fw, rv, window=window)
            data_by_name[name] = p_edited
    data_by_condition = {
        condition: {
            name: data_by_name[name][condition]
            for name in data_by_name
        }
        for condition in data_by_name[name]
    }
    return data_by_condition, run_names

data_path = sys.argv[1]
window = None
if len(sys.argv) > 2:
    window = int(sys.argv[2])
read_path = f"{data_path}/readfiles.csv"
condition_path = f"{data_path}/barcodes.csv"
sequence = read_sequence(f"{data_path}/sequence.fa")
guide = read_sequence(f"{data_path}/guide.fa")
out_path = f"{data_path}/outputs.csv"

condition_dict, index_list, condition_names = read_condition_dict(condition_path)
example_fw, example_rv = condition_dict[condition_names[0]]
barcode_length = len(example_fw)
data_by_condition, run_names = run_all_demux(condition_dict, index_list, condition_names, data_path, read_path,
                                             window=window)
with open(out_path, "wt") as f:
    f.write("condition name,index," + ",".join(run_names) + "\n")
    for index, condition in zip(index_list, condition_names):
        line = f"{condition},{index}," + ",".join([
            str(data_by_condition[condition][name]) for name in run_names])
        f.write(line + "\n")
