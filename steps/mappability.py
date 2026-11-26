#!/usr/bin/env python3
import argparse
import io
import subprocess
import sys
import threading
import re
import gzip
import bz2
from Bio import SeqIO
from intervaltree import IntervalTree, Interval
from collections import defaultdict

def open_file(filename, mode):
    if filename.endswith(".gz"):
        return gzip.open(filename, mode)
    elif filename.endswith(".bz2"):
        return bz2.open(filename, mode)
    else:
        return open(filename, mode)

def nonempty(value: str) -> str:
    if value == "":
        raise argparse.ArgumentTypeError("must not be empty")
    return value

def positive(value: str) -> int:
    try:
        number = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"invalid number '{value}'")
    if number <= 0:
        raise argparse.ArgumentTypeError(f"must be positive, got {number}")
    return number

def names(value: str) -> list:
    return [ name for name in value.split(",") if name != "" ]

# Command-line arguments
PARSER = argparse.ArgumentParser(
    description="Find the regions of a genome where reads cannot be mapped reliably, by simulating "
                "a read or read pair at every position, mapping them with bwa mem, and checking "
                "where they were mapped to. Arguments after -- are passed to bwa mem.")
PARSER.add_argument("--genome", type=nonempty, required=True, metavar="FILE",
                    help="genome sequences (FASTA)")
PARSER.add_argument("--genome-index", type=nonempty, required=True, metavar="PREFIX",
                    help="bwa index of the genome")
PARSER.add_argument("--circular", action="store_true",
                    help="sequences are circular, i.e. fragments may straddle the origin (default: linear)")
PARSER.add_argument("--paired", action="store_true",
                    help="reads are paired-end (default: single-end)")
PARSER.add_argument("--fragment-length", type=positive, metavar="N",
                    help="length of the simulated fragments (required with --paired, not allowed without)")
PARSER.add_argument("--read-length", type=positive, required=True, metavar="N",
                    help="length of the simulated reads")
PARSER.add_argument("--mapping-quality", type=int, required=True, metavar="MAPQ",
                    help="reads mapped with lower mapping quality count as not mappable")
PARSER.add_argument("--regions", type=nonempty, metavar="FILE",
                    help="regions (BED with names in the 4th column, optionally gzip- or bzip2-compressed) "
                         "that can be excluded with --exclude-regions")
PARSER.add_argument("--exclude-regions", type=names, default=[], metavar="NAMES",
                    help="comma-separated names of regions from --regions to mark as non-mappable")
PARSER.add_argument("--output", type=nonempty, required=True, metavar="FILE",
                    help="BED file of non-mappable regions to write")
PARSER.add_argument("--verbose", action="store_true",
                    help="report how each simulated read was mapped")
PARSER.add_argument("bwa_args", nargs="*", metavar="BWA_ARG",
                    help="additional arguments for bwa mem, after --")
ARGS = PARSER.parse_args()
if ARGS.paired and ARGS.fragment_length is None:
    PARSER.error("--fragment-length is required with --paired")
if not ARGS.paired and ARGS.fragment_length is not None:
    PARSER.error("--fragment-length is only allowed with --paired; for single-end reads, each simulated fragment is a read")
if ARGS.exclude_regions and ARGS.regions is None:
    PARSER.error("--exclude-regions requires --regions")

GENOME = ARGS.genome
GENOME_INDEX = ARGS.genome_index
CIRCULAR = ARGS.circular
READLENGTH = ARGS.read_length
IS_PAIRED = ARGS.paired
# For single-end reads, each simulated fragment is a single read
FRAGMENTLENGTH = ARGS.fragment_length if IS_PAIRED else READLENGTH
MAPQ = ARGS.mapping_quality
REGIONS = ARGS.regions or ""
EXCLUDE = ARGS.exclude_regions
DISPLACEMENT = 1
OUTPUT = ARGS.output
BWA_ARGS = ARGS.bwa_args
VERBOSE = ARGS.verbose

# Read genome to find sequence length and create interval tree
# containing mappable regions for each sequence
SEQLENGTHS = dict()
MAPABLE = dict()
for r in SeqIO.parse(GENOME, "fasta"):
    # Start with whole sequence marked as mappable
    seq = r.seq
    seqname = r.id
    seqlength = len(r.seq)
    SEQLENGTHS[seqname] = seqlength
    MAPABLE[seqname] = IntervalTree([Interval(0, seqlength)])

# Remove manually excluded regions for the interval tree
if REGIONS != "":
    bed_file = open_file(REGIONS, "rt")
    for line in bed_file:
        # Skip empty lines, comments and header lines
        if (line.strip() == "") or line.startswith(("#", "track", "browser")):
            continue
        fields = [ f.strip() for f in line.split("\t") ]
        if len(fields) < 4:
            raise RuntimeError(f"Regions BED file must have at least 4 columns, found only {len(fields)} columns in line: {line}")
        chr, start, end, name = fields[0], int(fields[1]), int(fields[2]), fields[3]
        # Resolve feature names
        if name in EXCLUDE:
            # Make sure min < max
            if start > end:
                start, end = end, start
            # Mark range as non-mappable
            if chr in MAPABLE:
                MAPABLE[chr].chop(start, end)
                print(f"Manually excluded region {name} at {chr}:{start+1}-{end}", file=sys.stderr)
            else:
                print(f"Manually excluded region {name} lies on unknown sequence {chr}", file=sys.stderr)

# Launch mapper
# For paired reads, tell bwa the insert size distribution instead of letting it
# infer one. All simulated fragments have the same length, so bwa would infer a
# standard deviation of zero, which prevents it from using the mate to place reads
# (e.g. within repeats). Inference can also make bwa accept other orientations
# (e.g. outward-facing RF pairs, as produced by fragments spanning the end of a
# circular sequence) as concordant; with -I, only FR pairs are.
# The standard deviation is set to 15% of the mean to mimic the width of a real
# insert size distribution, so that bwa is as uncertain about the placement of
# pairs as it would be for real reads (e.g. for nearby repeat copies). bwa accepts
# insert sizes within 4 standard deviations of the mean. Which pairs count as
# correctly mapped is decided below, using the exact TLEN.
# Note that bwa measures the insert size from the start of the forward read to the
# last base of the reverse read, i.e. it is the fragment length minus one (whereas
# the SAM TLEN equals the fragment length). BWA_ARGS come later, so they can
# override this.
BWA_INSERTSIZE = FRAGMENTLENGTH - 1
BWA_PAIRED_ARGS = ["-I", f"{BWA_INSERTSIZE},{0.15 * BWA_INSERTSIZE:.1f}"] if IS_PAIRED else []
BWA = subprocess.Popen(
    ["bwa", "mem"] + BWA_PAIRED_ARGS + BWA_ARGS + (["-p"] if IS_PAIRED else []) + [GENOME_INDEX, "-"],
    stdin=subprocess.PIPE, stdout=subprocess.PIPE,
    shell=False, encoding="ASCII"
)

# Read alignment produced by BWA and remove non-mappable regions
# from the per-chromosome interval trees we keep in MAPABLE. Reading
# alignment is done on a separate thread to avoid dead-lock issues.
def process_alignments():
    qname_pattern = re.compile(r"^(.*):([0-9]*):([0-9]*)$")
    is_prop, is_unmapped, is_rev, is_1st, is_2nd = 0x2, 0x4, 0x10, 0x40, 0x80
    for line in BWA.stdout:
        # Skip header lines
        if line[0] == "@":
            continue
        # Parse SAM line
        # Note: Convert from 1-based to 0-based coordinates because 1-based is confusing!
        qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, *tags = line.split("\t")
        flag, pos, mapq, pnext, tlen = int(flag), int(pos) - 1, int(mapq), int(pnext) - 1, int(tlen)
        # Extract true seq name, fragment start & end from query name
        m = qname_pattern.match(qname)
        if m is None:
            raise RuntimeError(f"Unexpected read name: {qname}")
        rname_true, fstart, fend = m[1], int(m[2]), int(m[3])
        # Determine whether the pair is concordant. bwa's proper-pair flag isn't
        # enough: bwa accepts a range of insert sizes, but all simulated fragments
        # have the same length, so a pair with any other insert size was mapped
        # wrongly. Only accept FR pairs with exactly the expected TLEN, i.e.
        # +FRAGMENTLENGTH for the forward and -FRAGMENTLENGTH for the reverse read.
        # This also prevents the range between the reads, which is flagged below as
        # the wrong target, from becoming arbitrarily large (e.g. for outward-facing
        # pairs, or for supplementary alignments, which inherit the flag but whose
        # TLEN is measured to the mate).
        is_concordant = (flag & is_prop != 0)
        if is_concordant and (tlen != (FRAGMENTLENGTH if flag & is_rev == 0 else -FRAGMENTLENGTH)):
            if VERBOSE:
                print(f"DISCORDANT: {qname} flagged as proper pair by bwa, but flag={flag} tlen={tlen}", file=sys.stderr)
            is_concordant = False
        # Check if the read was mapped correctly
        is_wrong = (rname != rname_true) or (mapq < MAPQ) or (cigar != f"{len(seq)}M")
        if IS_PAIRED:
            is_wrong = is_wrong or not is_concordant
            if (flag & is_1st != 0):
                is_wrong = is_wrong or (pos != fstart)
            elif (flag & is_2nd != 0):
                # The second read covers the last rlen bases of the fragment. Compute
                # rlen like the read generator does, i.e. reads are truncated if the
                # fragment is shorter than READLENGTH. If fend <= fstart, the fragment
                # straddles the origin of a circular sequence, and so may the read.
                seqlength_true = SEQLENGTHS[rname_true]
                flen = fend - fstart if fend > fstart else fend - fstart + seqlength_true
                rlen = min(READLENGTH, flen)
                is_wrong = is_wrong or (pos != (fend - rlen) % seqlength_true)
        else:
            is_wrong = is_wrong or (pos != fstart)
        if is_wrong:
            # Flag true fragment location as unmappable
            if VERBOSE:
                print(f"MISMATCH: {rname}:{fstart}-{fend} mapped at {rname}:{pos} flag={flag} ({cigar}, {mapq})", file=sys.stderr)
                print(f"        : WRONG SOURCE: {rname}:{fstart}-{fend}", file=sys.stderr)
            if fstart < fend:
                MAPABLE[rname_true].chop(fstart, fend)
            else:
                # Fragment straddles the origin of a circular sequence, i.e. it
                # covers [fstart, seqlength) and [0, fend)
                MAPABLE[rname_true].chop(fstart, SEQLENGTHS[rname_true])
                MAPABLE[rname_true].chop(0, fend)
            # Flag location to which the fragment was wrongly mapped as unmappable
            if IS_PAIRED and is_concordant:
                # For concordant pairs, use the range between the two reads.
                # |TLEN| is the number of bases from the leftmost mapped base
                # to the rightmost one, i.e. the (half-open) range has length |TLEN|.
                if tlen > 0:
                    # First read, range extends tlen bases from its start
                    if VERBOSE:
                        print(f"        : WRONG TARGET: {rname}:{pos}-{pos+tlen} (left mate)", file=sys.stderr)
                    MAPABLE[rname].chop(pos, pos + tlen)
                elif tlen < 0:
                    # Last read, range extends -tlen bases from its mate's start
                    if VERBOSE:
                        print(f"        : WRONG TARGET: {rname}:{pnext}-{pnext-tlen} (righ mate)", file=sys.stderr)
                    MAPABLE[rname].chop(pnext, pnext - tlen)
            else:
                # Otherwise, just use the range where the read mapped to
                # Note: Technically, we should compute the length along the reference from the CIGAR,
                # but we just use the sequence length instead for simplicity
                # Unmapped reads have no target to flag
                if flag & is_unmapped == 0:
                    if VERBOSE:
                        print(f"        : WRONG TARGET: {rname}:{pos}-{pos + len(seq)} (singleton)", file=sys.stderr)
                    MAPABLE[rname].chop(pos, pos + len(seq))
        else:
            if VERBOSE:
                print(f"   MATCH: {rname}:{fstart}-{fend} mapped as {rname}:{pos} flag={flag} ({cigar})", file=sys.stderr)
            pass

# Exceptions don't propagate from the reader thread to the main thread, so
# store them for the main thread to report. Also kill bwa, otherwise bwa blocks
# once its output isn't read anymore, and the main thread then blocks writing to
# bwa's input. Killing bwa makes these writes fail instead.
READER_ERROR = None
def read_alignments():
    global READER_ERROR
    try:
        process_alignments()
    except BaseException as e:
        READER_ERROR = e
        BWA.kill()
BWA_READER = threading.Thread(target=read_alignments, daemon=True)
BWA_READER.start()

# Close bwa's input, wait for bwa and the reader thread to finish, and fail
# if either of them failed. Without this check, a failure of bwa would yield
# incomplete results, e.g. an empty output (i.e. everything mappable).
def finish_bwa():
    try:
        BWA.stdin.close()
    except BrokenPipeError:
        pass
    BWA.wait()
    BWA_READER.join()
    if READER_ERROR is not None:
        raise RuntimeError("Failed to process the alignments produced by bwa") from READER_ERROR
    if BWA.returncode != 0:
        raise RuntimeError(f"bwa mem failed with exit code {BWA.returncode}")

# Send a read to bwa. If bwa stopped reading its input, find out why and fail.
def send_read(name, read):
    try:
        print(f">{name}\n{read}", file=BWA.stdin)
    except BrokenPipeError:
        finish_bwa()
        raise RuntimeError("bwa mem stopped reading its input unexpectedly")

# Input synthetic reads to BWA
for r in SeqIO.parse(GENOME, "fasta"):
    seq = r.seq
    seqname = r.id
    seqlength = len(r.seq)
    # Loop over fragment starting points. Starting points left of the sequence
    # start produce fragments truncated at the start of linear sequences.
    for fstart in range(-FRAGMENTLENGTH+1, seqlength, DISPLACEMENT):
        # Determine end point (exclusive)
        if CIRCULAR:
            # Let fragments straddle the origin. Starting points left of the
            # sequence start produce the same fragments as starting points at
            # the end of the sequence, so skip them.
            if (fstart < 0):
                continue
            # Wrap the (exclusive) end point into (0, seqlength], so a fragment
            # ending exactly at the genome end gets fend = seqlength (not 0).
            # For fragments straddling the origin, fend < fstart.
            fend = fstart + FRAGMENTLENGTH
            if fend > seqlength:
                fend -= seqlength
            flen = FRAGMENTLENGTH
        else:
            # Truncate fragments at sequence boundaries
            fend = min(fstart + FRAGMENTLENGTH, seqlength)
            fstart = max(fstart, 0)
            flen = fend - fstart
        # Generate reads
        # Ensure that the read length does not exceed the fragment length
        rlen = min(READLENGTH, flen)
        if fstart + rlen <= seqlength:
            # Read does not wrap around genome boundary
            r1 = str(seq[fstart:(fstart+rlen)].upper())
        else:
            # Read wraps around genome boundaries
            r1 = str(seq[fstart:seqlength].upper()) + \
                 str(seq[0:(fstart+rlen-seqlength)].upper())
        if not IS_PAIRED:
            # Single-end reads: the fragment is the read
            send_read(f"{seqname}:{fstart}:{fend}", r1)
            continue
        if fend - rlen >= 0:
            # Read does not wrap around genome boundary
            r2 = str(seq[(fend-rlen):fend].reverse_complement().upper())
        else:
            # Read wraps around genome boundaries
            r2 = str(seq[0:fend].reverse_complement().upper()) + \
                 str(seq[(fend-rlen+seqlength):seqlength].reverse_complement().upper())
        # Output read pair
        send_read(f"{seqname}:{fstart}:{fend}/1", r1)
        send_read(f"{seqname}:{fstart}:{fend}/2", r2)

# Done producing reads, close input, wait for BWA and reader to finish
finish_bwa()

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


