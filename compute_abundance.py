#!/usr/bin/env python3
import yaml
import argparse
import sys
import os

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

print("*** MAP {sample} -> {bam_file}".format_map(CFG), file=sys.stderr)
if not os.path.exists(CFG["bam_file"]):
    print(CFG["map_cmd"])
else:
    print(f"  {bam_file} already exists, reusing...")

