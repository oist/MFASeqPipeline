#!/usr/bin/env python3
import io
import subprocess
import sys
import threading
import re
from Bio import SeqIO
from intervaltree import IntervalTree, Interval
from collections import defaultdict

# Command-line arguments
GENOME = sys.argv[1]
CIRCULAR = bool(sys.argv[2])
FRAGMENTLENGTH = int(sys.argv[3])
READLENGTH = int(sys.argv[4])
DISPLACEMENT = 1
OUTPUT = sys.argv[5]
BWA_ARGS = sys.argv[6:]
VERBOSE = False

BWA = subprocess.Popen(
    ["bwa", "mem"] + BWA_ARGS + ["-p"] + [GENOME, "-"],
    stdin=subprocess.PIPE, stdout=subprocess.PIPE,
    shell=False, encoding="ASCII"
)

# Read alignment produced by BWA and remove non-mappable reagions
# from the per-chromosome interval trees we keep in MAPABLE. Reading
# alignment is done on a separate thread to avoid dead-lock issues.
MAPABLE = dict()
def read_alignments():
    qname_pattern = re.compile(r"^(.*):([0-9]*):([0-9]*)$")
    is_prop, is_rev, is_1st, is_2nd = 0x2, 0x10, 0x40, 0x80
    for line in BWA.stdout:
        # Skip header lines
        if line[0] == "@":
            continue
        # Parse SAM line
        # Note: Convert from 1-based to 0-based coordinates because 1-based is confusing!
        qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, *tags = line.split("\t")
        flag, pos, pnext, tlen = int(flag), int(pos) - 1, int(pnext) - 1, int(tlen)
        # Extract true seq name, fragment start & end from query name
        m = qname_pattern.match(qname)
        if m is None:
            raise RuntimeError(f"Unexpected read name: {qname}")
        rname_true, fstart, fend = m[1], int(m[2]), int(m[3])
        # Check if the read was mapped correctly
        is_wrong = (flag & is_prop == 0) or (rname != rname_true) or (cigar != f"{len(seq)}M")
        if (flag & is_1st != 0):
            is_wrong = is_wrong or (pos != fstart)
        elif (flag & is_2nd != 0):
            is_wrong = is_wrong or (pos != fend - READLENGTH)
        if is_wrong:
            # Flag true fragment location as unmappable
            if VERBOSE:
                print(f"MISMATCH: {rname}:{fstart}-{fend} mapped at {rname}:{pos} flag={flag} ({cigar})", file=sys.stderr)
                print(f"        : WRONG SOURCE: {rname}:{fstart}-{fend}", file=sys.stderr)
            MAPABLE[rname_true].chop(fstart, fend)
            # Flag location to which the fragment was wrongly mapped as unmappable
            if (flag & is_prop != 0):
                # For proper pairs, use the range between the two reads
                if tlen > 0:
                    # First read, range extends tlen further bases
                    if VERBOSE:
                        print(f"        : WRONG TARGET: {rname}:{pos}-{pos+tlen+1} (left mate)", file=sys.stderr)
                    MAPABLE[rname_true].chop(pos, pos + tlen + 1)
                elif tlen < 0:
                    # Last read, range extends tlen bases from first
                    if VERBOSE:
                        print(f"        : WRONG TARGET: {rname}:{pnext}-{pnext-tlen+1} (righ mate)", file=sys.stderr)
                    MAPABLE[rname_true].chop(pnext, pnext + (-tlen) + 1)
            else:
                # Otherwise, just use the range where the read mapped to
                # Note: Technically, we should compute the length along the reference from the CIGAR,
                # but we just use the sequene length instead for simplicity
                if VERBOSE:
                    print(f"        : WRONG TARGET: {rname}:{pos}-{pos + len(seq)} (singleton)", file=sys.stderr)
                    MAPABLE[rname_true].chop(pos, pos + len(seq))
        else:
            if VERBOSE:
                print(f"   MATCH: {rname}:{fstart}-{fend} mapped as {rname}:{pos} flag={flag} ({cigar})", file=sys.stderr)
            pass
BWA_READER = threading.Thread(target=read_alignments, daemon=True)
BWA_READER.start()

# Input synthetic reads to BWA
SEQLENGTHS = dict()
for r in SeqIO.parse(GENOME, "fasta"):
    # Start with whole sequence marked as mappable
    seq = r.seq
    seqname = r.id
    seqlength = len(r.seq)
    SEQLENGTHS[seqname] = seqlength
    MAPABLE[seqname] = IntervalTree([Interval(0, seqlength)])
    # Loop over fragment starting points
    for fstart in range(-FRAGMENTLENGTH+1, seqlength-1, DISPLACEMENT):
        # Mark whole sequence as mapable
        # Determine end point (exclusive)
        if CIRCULAR:
            # Let fragments wrap around genome boundaries. Only
            # necessary for the end point, though. Starting points left
            # of genome start produce the same fragments as starting points
            # at the very end of the genome.
            # fragments as starting points at the very end in this case.
            fend = (fstart + FRAGMENTLENGTH) % seqlength
            if (fstart < 0):
                continue
        else:
            # Truncate fragments at genome boundaries
            fend = min(fstart + FRAGMENTLENGTH, seqlength)
            fstart = max(fstart, 0)
        # Generate reads
        # Ensure that the read length does not exceed the fragment length
        flen = fend - fstart
        rlen = min(READLENGTH, flen)
        if fstart + rlen <= seqlength:
            # Read does not wrap around geome boundary
            r1 = str(seq[fstart:(fstart+rlen)].upper())
        else:
            # Read wraps around genome boundaries
            r1 = str(seq[fstart:seqlength].upper()) + \
                 str(seq[0:(fstart+rlen-seqlength)].upper())
        if fend - rlen >= 0:
            # Read does not wrap around genome boundary
            r2 = str(seq[(fend-rlen):fend].reverse_complement().upper())
        else:
            # Read wraps around genome boundaries
            r2 = str(seq[0:fend].reverse_complement().upper()) + \
                 str(seq[(fend-rlen+seqlength):seqlength].reverse_complement().upper())
        # Output read pair                 
        print(f">{seqname}:{fstart}:{fend}/1\n{r1}", file=BWA.stdin)
        print(f">{seqname}:{fstart}:{fend}/2\n{r2}", file=BWA.stdin)

# Done producing reads, close input, wait for BWA and reader to finish
BWA.stdin.close()
BWA.wait()
BWA_READER.join()

# Output non-mappable regions
with open(OUTPUT, "w+", encoding="ASCII") as out:
    print(f"# Command used to create this file:", file=out)
    print(f"#   {' '.join(sys.argv)}", file=out)
    print(f"# Chromosomes:", file=out)
    for seqname, length in SEQLENGTHS.items():
        print(f"#   {seqname}: {length}", file=out)
    # Iterate over sequences
    for seqname, mappable in MAPABLE.items():
        # Create interval tree of non-mappable regions from tree of mappable regions
        nonmappable = IntervalTree([Interval(0, SEQLENGTHS[seqname])])
        for range in mappable:
            nonmappable.chop(range.begin, range.end)
        # Output non-mappable regions
        for range in sorted(nonmappable):
            print(f"{seqname}\t{range.begin}\t{range.end}", file=out)


