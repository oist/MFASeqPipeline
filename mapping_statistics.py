#!/usr/bin/env python3
import sys
import statistics
import pysam
from collections import namedtuple

bam_file = sys.argv[1]
stats_file = sys.argv[2]

Stats = namedtuple("Stats", ["total_reads", "mapped_reads", "properly_paired_reads", "insert_size", "read_length"])

bam = pysam.AlignmentFile(bam_file, mode="rb")
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
with open(stats_file, "w+") as f:
    for field, value in stats._asdict().items():
        print(f"{field}\t{value}", file=f)
