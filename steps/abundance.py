#!/usr/bin/env python3
import sys
import os
import io
import gzip
import math
import numpy as np
import pysam
from intervaltree import IntervalTree, Interval

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "ext", "pyngstuff"))
import pyngstuff

class Coverage:
    def __init__(self, length : int):
        self.length = length
        self.deltas = np.zeros([length], dtype=np.int32)
        self.startingpoints = np.zeros([length], dtype=np.int32)
        self.endpoints = np.zeros([length], dtype=np.int32)
        self.fragment_length_sum = 0
        self.fragment_length_count = 0

    def fragment_inclusive(self, start : int, end : int):
        assert (start >= 0) and (start < self.length)
        assert (end >= 0) and (end < self.length)
        assert start <= end
        # Mark fragment start- and endpoint
        self.startingpoints[start] += 1
        self.endpoints[end] += 1
        # Increase coverage on [start, end]
        self.deltas[start] += 1
        if end < self.length - 1:
            self.deltas[end+1] -= 1
        self.fragment_length_sum += (end - start + 1)
        self.fragment_length_count += 1

    def coverage(self) -> np.ndarray:
        return np.cumsum(self.deltas)
        
    def fragment_length(self) -> float:
        if self.fragment_length_count > 0:
            return round(self.fragment_length_sum / self.fragment_length_count)
        else:
            return math.nan
        
def fragment_coverage(wsize: int, flen_float: float) -> float:
    # Compute coverage length, i.e. average number of bases of a fragment that covers a window
    # For fragments of length k and the window [0, n), fragments that overlap the window start at
    # a base in (-k, n) = [-k+1, n-1]. A fragment that starts at base i covers the bases
    # [i, i + k - 1], and overlaps the window on [ max(i, 0), min(i + k - 1, n - 1) ]. The
    # expected coverage (assuming uniformly distributed starting points) is thus
    #      [ sum i=-k+1 to n-1 min(i + k, n) - max(i, 0) ] / (n - k + 1)
    if not math.isfinite(flen_float):
        return math.nan
    flen: int = int(round(flen_float))
    return sum([ min([i + flen, wsize]) - max([i, 0])
                 for i in range(-flen+1, wsize)]) / (wsize + flen - 1)

def window_statistics(coverage: Coverage, wsize: int) -> np.ndarray:
    # Determine padding to append to make length divisible by window size
    cov = coverage.coverage()
    pad = (-cov.shape[0]) % wsize
    # Reshape coverage vector into a matrix where each column is a window
    cov_pad = np.concatenate((cov, np.zeros(pad, dtype=cov.dtype)))
    cov_win = np.reshape(cov_pad, (wsize, -1), order='F')
    # Reshape startingpoint vector into a matrix where each column is a window
    sp_pad = np.concatenate((coverage.startingpoints, np.zeros(pad, dtype=cov.dtype)))
    sp_win = np.reshape(sp_pad, (wsize, -1), order='F')
    ep_pad = np.concatenate((coverage.endpoints, np.zeros(pad, dtype=cov.dtype)))
    ep_win = np.reshape(ep_pad, (wsize, -1), order='F')
    # Create matrix of positions with the same shape
    i_win = np.reshape(np.arange(start=0, stop=cov_win.size), (wsize, -1), order='F')
    # Compute statistics for each window
    first = np.min(i_win, axis=0)
    last = np.max(i_win, axis=0)
    mid = np.average(i_win, axis=0)
    wlen = np.full(cov_win.shape[1], fill_value=wsize)
    flen = np.full(cov_win.shape[1], fill_value=coverage.fragment_length())
    fcov = np.full(cov_win.shape[1], fill_value=fragment_coverage(wsize, flen[0]))
    cov_avg = np.average(cov_win, axis=0)
    cov_sd = np.std(cov_win, axis=0)
    sp = np.sum(sp_win, axis=0)
    ep = np.sum(ep_win, axis=0)
    # Fix last window
    if pad > 0:
        i_lastwindow = i_win[0:-pad, -1]
        cov_lastwindow = cov_win[0:-pad, -1]
        first[-1] = np.min(i_lastwindow)
        last[-1] = np.max(i_lastwindow)
        mid[-1] = np.average(i_lastwindow)
        wlen[-1] = cov_lastwindow.shape[0]
        fcov[-1] = fragment_coverage(wlen[-1], flen[-1])
        cov_avg[-1] = np.average(cov_lastwindow)
        cov_sd[-1] = np.std(cov_lastwindow)
    # We estimate the standard error of the window average by interpreting it as a
    # transformed Poissonian random variable. The average coverage A (= avg) is
    # is approximately N * L / W where N is the number of fragments overlapping the
    # window, L (= covlen) the average number of bases covered by a single fragment,
    # and W (= wlen) the window size. N should approximately follow a Poissonian
    # distribution and thus have standard error sqrt(N). Therefore the standard error
    # of A = N * L / W is
    #   sqrt(N) * L / W = sqrt(A * W / L) * L / W = sqrt(A) * sqrt(L/W)
    cov_serr = np.sqrt(cov_avg) * np.sqrt(fcov / wlen)
    # Create result table
    cols = [('first', np.int32), ('mid', np.float64), ('last', np.int32),
            ('wlen', np.int32), ('flen', np.float64), ('fcov', np.float64),
            ('avg', np.float64), ('sd', np.float64), ('serr', np.float64), 
            ('sp', np.int32), ('ep', np.int32)]
    tab = np.zeros(wlen.shape[0], dtype=cols)
    tab["first"] = first
    tab["mid"] = mid
    tab["last"] = last
    tab["wlen"] = wlen
    tab["flen"] = flen
    tab["fcov"] = fcov
    tab["avg"] = cov_avg
    tab["sd"] = cov_sd
    tab["serr"] = cov_serr
    tab["sp"] = sp
    tab["ep"] = ep
    return tab

def ratio(sample: np.ndarray, ref: np.ndarray) -> np.ndarray:
    # To estimate the standard error of ratios, we use
    #
    #                     CV[X]^2 + CV[Y]^2 - CV[X]^2 CV[Y]^2
    #    CV[ X / Y ]^2 = -------------------------------------
    #                               1 - 2 * CV[Y]^2 
    #
    # where CV[X] = STDDEV[X] / MEAN[X] denotes the coefficient of variation
    # 
    # This identity holds exactly only for Gamma-distributed random
    # variables, but should hold approximately for any RV with short tails
    # We compute the ratio of average coverage, startingpoints and endpoints
    # between sample and reference
    cols = [('avg', np.float64), ('avg_serr', np.float64),
            ('sp', np.float64), ('sp_serr', np.float64),
            ('ep', np.float64), ('ep_serr', np.float64)]
    tab = np.zeros(sample.shape[0], dtype=cols)
    # Ratio of average coverages
    tab["avg"] = sample["avg"] / ref["avg"]
    s_cv2 = (sample["serr"] / sample["avg"]) ** 2
    r_cv2 = (ref["serr"] / ref["avg"]) ** 2
    q_cv2 = (s_cv2 + r_cv2 - s_cv2 * r_cv2) / (1 - 2 * r_cv2)
    tab["avg_serr"] = tab["avg"] * np.sqrt(q_cv2)
    # Ratio of startingpoints
    tab["sp"] = sample["sp"] / ref["sp"]
    s_cv2 = 1 / sample["sp"]
    r_cv2 = 1 / ref["sp"]
    q_cv2 = (s_cv2 + r_cv2 - s_cv2 * r_cv2) / (1 - 2 * r_cv2)
    tab["sp_serr"] = tab["sp"] * np.sqrt(q_cv2)
    # Ratio of endpoints
    tab["ep"] = sample["ep"] / ref["ep"]
    s_cv2 = 1 / sample["ep"]
    r_cv2 = 1 / ref["ep"]
    q_cv2 = (s_cv2 + r_cv2 - s_cv2 * r_cv2) / (1 - 2 * r_cv2)
    tab["ep_serr"] = tab["ep"] * np.sqrt(q_cv2)
    return tab
    

def read_bam(bamfile, is_paired, mapq):
    with pysam.AlignmentFile(bamfile, "rb") as bam:
        # Read sequence lengths from BAM header
        seqlengths = { seq: length for seq, length in zip(bam.header.references, bam.header.lengths) }
        # Create Coverage objects for all chromosomes
        cov = {seq: Coverage(length) for seq, length in seqlengths.items()}
        # Scan read pairs
        if is_paired:
            # Scan read pairs, mark whole fragment as covered
            pairs = pyngstuff.read_pair_iterator(bam)
            for pair in pairs:
                if (pair.forward.mapping_quality < mapq) or (pair.reverse.mapping_quality < mapq):
                    continue
                if pair.forward.reference_start >= pair.reverse.reference_end:
                    continue
                cov[pair.forward.reference_name].fragment_inclusive(
                    start=pair.forward.reference_start, end=pair.reverse.reference_end - 1)
        else:
            # Scan single reads, mark only the read as covered
            for read in bam:
                if read.is_unmapped:
                    continue
                if read.is_secondary or read.is_supplementary:
                    continue
                if read.mapping_quality < mapq:
                    continue
                cov[read.reference_name].fragment_inclusive(
                  start=read.reference_start, end=read.reference_end - 1)

        return (cov, seqlengths)


def read_mappability(mappabilityfile, seqlength):
    mappability = { seq: IntervalTree() for seq, length in seqlength.items() }
    with open(mappabilityfile, "r") as f:
        for line in f.readlines():
            # Skip comments
            if line.lstrip().startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                raise RuntimeError(f"too few fields in line {line}")
            seq, start, end = fields[0], int(fields[1]), int(fields[2])
            mappability[seq].add(Interval(start, end))
    return mappability

def open_output(filename: str):
    if filename.endswith(".gz"):
        return gzip.open(filename, "wt+")
    else:
        return open(filename, "wt+")

# Make output unbuffered
sys.stdout = io.TextIOWrapper(open(sys.stdout.fileno(), 'wb', 0), write_through=True)

# Command-line arguments
SAM_BAM_FILE = sys.argv[1]
SAM_MAPPABILITY_FILE = sys.argv[2]
REF_BAM_FILE = sys.argv[3]
REF_MAPPABILITY_FILE = sys.argv[4]
IS_PAIRED = sys.argv[5] in ["1", "true"]
MAPQ = int(sys.argv[6])
WINDOWSIZE = sys.argv[7]
ABUNDANCE_FILE = sys.argv[8]

if WINDOWSIZE.endswith("k"):
    WINDOWSIZE = int(WINDOWSIZE[0:-1]) * 1000
else:
    WINDOWSIZE = int(WINDOWSIZE)

# Read data
PAIRED_DESC = "paired-end" if IS_PAIRED else "single-end"
print(f"Reading {PAIRED_DESC} sample BAM {SAM_BAM_FILE}, using reads with mapping quality >= {MAPQ}")
SAM_COVERAGE, SAM_SEQLENGTHS = read_bam(SAM_BAM_FILE, IS_PAIRED, MAPQ)

if REF_BAM_FILE != "":
    print(f"Reading {PAIRED_DESC} reference BAM {REF_BAM_FILE}, using reads with mapping quality >= {MAPQ}")
    REF_COVERAGE, REF_SEQLENGTHS = read_bam(REF_BAM_FILE, IS_PAIRED, MAPQ)
else:
    print(f"No reference BAM file provided")
    REF_COVERAGE, REF_SEQLENGTHS = None, None

print(f"Reading {SAM_MAPPABILITY_FILE}")
SAM_MAPPABILITY = read_mappability(SAM_MAPPABILITY_FILE, SAM_SEQLENGTHS)
if REF_MAPPABILITY_FILE != "":
    print(f"Reading {REF_MAPPABILITY_FILE}")
    REF_MAPPABILITY = read_mappability(REF_MAPPABILITY_FILE, REF_SEQLENGTHS)
else:
    print(f"No reference mappability file provided")
    REF_MAPPABILITY = None


# Compute
#
# NOTE: In the following, the data is kept in dictionaries which contain
# one numpy.ndarray per chromosome, indexed by the chromosome name. 


# Compute average coverage and std. dev. within each window
# window_average() returns a heterogeonous ndarray with named columns
print(f"Computing abundances and abundance ratios")
SAM = { seq: window_statistics(SAM_COVERAGE[seq], WINDOWSIZE)
        for seq in SAM_SEQLENGTHS.keys() }
if REF_COVERAGE != None and REF_SEQLENGTHS != None:
    REF = { seq: window_statistics(REF_COVERAGE[seq], WINDOWSIZE)
            for seq in REF_SEQLENGTHS.keys() }
else:
    REF = None

# Compute coverage ratio
if REF != None:
    RATIO = { seq: ratio(SAM[seq], REF[seq]) for seq in SAM.keys() }
else:
    RATIO = None

# Test whether windows overlap some non-mappable regions
print(f"Finding mappable windows for window size {WINDOWSIZE}bp")
SAM_MAPPABLE = { seq: np.array([ not SAM_MAPPABILITY[seq].overlaps(start, end+1)
                                 for (start, end) in zip(win["first"].tolist(), win["last"].tolist()) ])
                 for seq, win in SAM.items() }
if REF != None and REF_MAPPABILITY != None:
    REF_MAPPABLE = { seq: np.array([ not REF_MAPPABILITY[seq].overlaps(start, end+1)
                                     for (start, end) in zip(win["first"].tolist(), win["last"].tolist()) ])
                     for seq, win in REF.items() }
else:
    REF_MAPPABLE = None

# Output results
print(f"Writing results to {ABUNDANCE_FILE}")
SLEN = max([len(seq) for seq in SAM_SEQLENGTHS.keys() ])
COLS = [('seq', f"U{SLEN}"), ('first', np.int32), ('midpoint', np.float64), ('last', np.int32),
        ('windowsize', np.int32), ('fragmentlength', np.float64), ('fragmentcoverage', np.float64), 
        ("sample_coverage_avg", np.float64), ("sample_coverage_sd", np.float64), ("sample_coverage_avg_stderr", np.float64),
        ("sample_startingpoints", np.int32), ("sample_endpoints", np.int32), ("sample_mappable", np.bool_), 
        ("reference_coverage_avg", np.float64), ("reference_coverage_sd", np.float64), ("reference_coverage_avg_stderr", np.float64),
        ("reference_startingpoints", np.int32), ("reference_endpoints", np.int32), ("reference_mappable", np.bool_),
        ("coverage_ratio", np.float64), ("coverage_ratio_stderr", np.float64),
        ("startingpoint_ratio", np.float64), ("startingpoint_ratio_stderr", np.float64),
        ("endpoint_ratio", np.float64), ("endpoint_ratio_stderr", np.float64),
        ("ratio_mappable", np.bool_) ]
FMT = ("%s",) + ("%d", "%.1f", "%d", "%d") + ("%.2f",)*2 + ("%.2f", "%.5f", "%.5f", "%d", "%d", "%d")*2 + ("%.5f", "%.5f")*3 + ("%d",)

with open_output(ABUNDANCE_FILE) as f:
    f.write("\t".join([name for (name, type) in COLS]) + "\n")
    for seq in SAM.keys():
        tab = np.zeros(SAM[seq].shape[0], dtype=COLS)
        tab["seq"] = np.full(SAM[seq].shape[0], fill_value=seq, dtype=tab["seq"].dtype)
        tab["first"] = SAM[seq]["first"]
        tab["midpoint"] = SAM[seq]["mid"]
        tab["last"] = SAM[seq]["last"]
        tab["windowsize"] = SAM[seq]["wlen"]
        tab["fragmentlength"] = SAM[seq]["flen"]
        tab["fragmentcoverage"] = SAM[seq]["fcov"]
        tab["sample_mappable"] = SAM_MAPPABLE[seq]
        tab["sample_coverage_avg"] = SAM[seq]["avg"]
        tab["sample_coverage_sd"] = SAM[seq]["sd"]
        tab["sample_coverage_avg_stderr"] = SAM[seq]["serr"]
        tab["sample_startingpoints"] = SAM[seq]["sp"]
        tab["sample_endpoints"] = SAM[seq]["ep"]
        if REF != None:
            tab["reference_coverage_avg"] = REF[seq]["avg"]
            tab["reference_coverage_sd"] = REF[seq]["sd"]
            tab["reference_coverage_avg_stderr"] = REF[seq]["serr"]
            tab["reference_startingpoints"] = REF[seq]["sp"]
            tab["reference_endpoints"] = REF[seq]["ep"]
        else:
            tab["reference_coverage_avg"] = np.nan
            tab["reference_coverage_sd"] = np.nan
            tab["reference_coverage_avg_stderr"] = np.nan
            tab["reference_startingpoints"] = -1
            tab["reference_endpoints"] = -1
        if REF_MAPPABLE != None:
            tab["reference_mappable"] = REF_MAPPABLE[seq]
        else:
            tab["reference_mappable"] = False
        if RATIO != None:
            tab["coverage_ratio"] = RATIO[seq]["avg"]
            tab["coverage_ratio_stderr"] = RATIO[seq]["avg_serr"]
            tab["startingpoint_ratio"] = RATIO[seq]["sp"]
            tab["startingpoint_ratio_stderr"] = RATIO[seq]["sp_serr"]
            tab["endpoint_ratio"] = RATIO[seq]["ep"]
            tab["endpoint_ratio_stderr"] = RATIO[seq]["ep_serr"]
        else:
            tab["coverage_ratio"] = np.nan
            tab["coverage_ratio_stderr"] = np.nan
            tab["startingpoint_ratio"] = np.nan
            tab["startingpoint_ratio_stderr"] = np.nan
            tab["endpoint_ratio"] = np.nan
            tab["endpoint_ratio_stderr"] = np.nan
        if REF_MAPPABLE != None:
            tab["ratio_mappable"] = SAM_MAPPABLE[seq] & REF_MAPPABLE[seq]
        else:
            tab["ratio_mappable"] = False
        np.savetxt(f, tab, delimiter="\t", fmt=FMT)
