#!/usr/bin/env python3
import yaml
import argparse
import sys
import os
import subprocess
import statistics
import pysam
from collections import namedtuple
#import pyngstuff

parser = argparse.ArgumentParser(
    prog="MFASeqPipeline", 
    description="Process MFA-Seq data")

parser.add_argument("config", help="Config file in YAML format")
parser.add_argument("sample", help="Sample name")
parser.add_argument("--fastq_dir", help="Directory containing FASTQ files")
parser.add_argument("--out_dir", help="Directory containing generated files")
parser.add_argument("--tmp_dir", help="Directory containing temporary files")
parser.add_argument("--genome", help="Reference genome")
parser.add_argument("--cores", help="Number of cores to use")
parser.add_argument("--memory", help="Memory to use (in gigabytes)")
parser.add_argument("--fragmentlength", help="Average fragment length (in bases)")
parser.add_argument("--resolution", help="Resolution of abundance data output (in bases)")

Stats = namedtuple("Stats", ["total_reads", "mapped_reads", "properly_paired_reads", "insert_size", "read_length"])

# Parse command-line arguments
ARGS = parser.parse_args()
with open(ARGS.config, "r") as cfg:
    CFG = yaml.safe_load(cfg)

# Command-line arguments override config
for key, value in vars(ARGS).items():
    if value is not None:
        CFG[key] = value

# Recursively expand {...} in config
for key in CFG.keys():
    v = CFG[key]
    if isinstance(v, str):
        vp = v.format_map(CFG)
        while v != vp:
            v = vp
            vp = v.format_map(CFG)
        CFG[key] = vp

# *** MAPPING
print("*** MAP {sample} -> {bam_file}".format_map(CFG), file=sys.stderr)
if not os.path.exists(CFG["bam_file"]):
    print("   Genome: {genome_index}".format_map(CFG))
    print("   Inputs:".format(CFG))
    print("     1st read: {read1_file}".format_map(CFG))
    print("     2nd read: {read2_file}".format_map(CFG))
    print("   Outputs:".format_map(CFG))
    print("     BAM: {bam_file}".format_map(CFG))
    subprocess.run(CFG["map_cmd"], shell=True, executable=CFG["shell"], stderr=sys.stdout, check=True)
else:
    print("  {bam_file} already exists, reusing...".format_map(CFG))

# *** INDEXING
print("*** INDEX {bam_file} -> {bam_index_file}".format_map(CFG), file=sys.stderr)
if not os.path.exists(CFG["bam_index_file"]):
    print("   Inputs:".format_map(CFG))
    print("     BAM: {bam_file}".format_map(CFG))
    print("   Outputs:".format_map(CFG))
    print("     index: {bam_index_file}".format_map(CFG))
    subprocess.run(CFG["index_cmd"], shell=True, executable=CFG["shell"], stderr=sys.stdout, check=True)
else:
    print("  {bam_index_file} already exists, reusing...".format_map(CFG))

# *** MAPPING STATISTICS
print("*** STATS {bam_file} -> {stats_file}".format_map(CFG), file=sys.stderr)
if not os.path.exists(CFG["stats_file"]):
    print("   Inputs:".format_map(CFG))
    print("     BAM: {bam_file}".format_map(CFG))
    print("   Outputs:".format_map(CFG))
    print("     Stats: {stats_file}".format_map(CFG))

    # Compute statistics
    bam = pysam.AlignmentFile(CFG["bam_file"], mode="rb")
    total_reads = 0
    mapped_reads = 0
    properly_paired_reads = 0
    read_lengths = []
    insert_sizes = []
    for read in bam.fetch():
        total_reads += 1
        read_lengths.append(read.query_length)
        if read.is_mapped:
            mapped_reads += 1
        if read.is_proper_pair:
            properly_paired_reads += 1
            insert_sizes.append(abs(read.isize))
    stats = Stats(total_reads, mapped_reads, properly_paired_reads,
                  int(statistics.median(insert_sizes)),
                  int(statistics.median(read_lengths)))

    # Write results to file
    with open(CFG["stats_file"], "w+") as f:
        for field, value in stats._asdict().items():
            print(f"{field}\t{value}", file=f)
else:
    print("  {stats_file} already exists, reusing...".format_map(CFG))
    stats_dict = dict()
    with open(CFG["stats_file"], "r") as f:
        for line in f.readlines():
            fields = line.split("\t", maxsplit=2)
            stats_dict[fields[0]] = int(fields[1])
    stats = Stats(**stats_dict)

# *** MAPPABILITY ANALYSIS
print("*** MAPPABILITY {bam_file} -> {mappability_file}".format_map(CFG), file=sys.stderr)
if not os.path.exists(CFG["mappability_file"]):
    print("   Inputs:".format_map(CFG))
    print("     BAM: {bam_file}".format_map(CFG))
    print("   Outputs:".format_map(CFG))
    print("     mappability: {mappability_file}".format_map(CFG))
    subprocess.run(CFG["mappability_cmd"], shell=True, executable=CFG["shell"], stderr=sys.stdout, check=True)
else:
    print("  {mappability_file} already exists, reusing...".format_map(CFG))

# *** ABUNDANCE COMPUTATION
print("*** ABUNDANCE {bam_file}, {mappability_file} -> {abundance_file}".format_map(CFG), file=sys.stderr)
if not os.path.exists(CFG["abundance_file"]):
    print("   Inputs:".format_map(CFG))
    print("     BAM: {bam_file}".format_map(CFG))
    print("     mappability: {mappability_file}".format_map(CFG))
    print("   Outputs:".format_map(CFG))
    print("     abundance: {abundance_file}".format_map(CFG))
    subprocess.run(CFG["abundance_cmd"], shell=True, executable=CFG["shell"], stderr=sys.stdout, check=True)
else:
    print("  {mappability_file} already exists, reusing...".format_map(CFG))
