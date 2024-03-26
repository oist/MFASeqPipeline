#!/usr/bin/env python3
import sys
import io
import numpy as np
import pyngstuff
import pysam
from intervaltree import IntervalTree, Interval

class Coverage:
    def __init__(self, length : int):
        self.length = length
        self.deltas = np.zeros([length], dtype=np.int32)
        self.fragment_length_sum = 0
        self.fragment_length_count = 0

    def fragment_inclusive(self, begin : int, end : int, dir : int):
        assert (begin >= 0) and (begin < self.length)
        assert (end >= 0) and (end < self.length)
        assert dir in [-1, 1]
        if dir == -1:
            return self.fragment_inclusive(end, begin, +1)
        if begin <= end:
            # [begin..end]
            self.deltas[begin] += 1
            if end < self.length - 1:
                self.deltas[end+1] -= 1
            self.fragment_length_sum += (end - begin + 1)
            self.fragment_length_count += 1
        elif (begin > end) and (dir == 1):
            # [0..end] u [begin...last]
            self.deltas[0] += 1
            self.deltas[end+1] -= -1
            self.deltas[begin] += 1
            self.fragment_length_sum += self.length - (begin - end - 1)
            self.fragment_length_count += 1

    def coverage(self) -> np.ndarray:
        return np.cumsum(self.deltas)
        
    def fragment_length(self) -> float:
        return int(self.fragment_length_sum / self.fragment_length_count)
        
def fragment_coverage(wsize: int, flen: int) -> float:
    # Compute coverage length, i.e. average number of bases of a fragment that covers a window
    # For fragments of length k and the window [0, n), fragments that overlap the window start at
    # a base in (-k, n) = [-k+1, n-1]. A fragment that starts at base i covers the bases
    # [i, i + k - 1], and overlaps the window on [ max(i, 0), min(i + k - 1, n - 1) ]. The
    # expected coverage (assuming uniformly distributed starting points) is thus
    #      [ sum i=-k+1 to n-1 min(i + k, n) - max(i, 0) ] / (n - k + 1)
    return sum([ min([i + flen, wsize]) - max([i, 0])
                 for i in range(-flen+1, wsize)]) / (wsize + flen - 1)

def window_statistics(cov: Coverage, wsize: int) -> np.ndarray:
    # Get coverage and fragment length
    v = cov.coverage()
    assert v.ndim == 1
    # Pad to make the number of entries a multiple of wsize
    pad = (-v.shape[0]) % wsize
    v_pad = np.concatenate((v, np.zeros(pad, dtype=v.dtype)))
    # Reshape into a matrix where each column is a window
    v_win = np.reshape(v_pad, (wsize, -1), order='F')
    # Compute statistics for each window
    wlen = np.full(v_win.shape[1], fill_value=wsize)
    flen = np.full(v_win.shape[1], fill_value=cov.fragment_length())
    fcov = np.full(v_win.shape[1], fill_value=fragment_coverage(wsize, flen[0]))
    avg = np.average(v_win, axis=0)
    sd = np.std(v_win, axis=0)
    # Fix last window
    if pad > 0:
        last = v_win[0:-pad, -1]
        wlen[-1] = last.shape[0]
        fcov[-1] = fragment_coverage(wlen[-1], flen[-1])
        avg[-1] = np.average(last)
        sd[-1] = np.std(last)
    # We estimate the standard error of the window average by interpreting it as a
    # transformed Poissonian random variable. The average coverage A (= avg)
    # is thenumber of fragments N overlapping the window times the average number of
    # bases L (= covlen) covered by a single fragment. We assume that the standard error
    # of N is sqrt(N), therefore the standard error of A = N * L is
    #   sqrt(N) * L = sqrt(A / L) * L = sqrt(A) * sqrt(L)
    serr = np.sqrt(avg) * np.sqrt(fcov)
    # Create result table
    cols = [('wlen', np.int32), ('flen', np.float64), ('fcov', np.float64),
            ('avg', np.float64), ('sd', np.float64), ('serr', np.float64) ]
    tab = np.zeros(wlen.shape[0], dtype=cols)
    tab["wlen"] = wlen
    tab["flen"] = flen
    tab["fcov"] = fcov
    tab["avg"] = avg
    tab["sd"] = sd
    tab["serr"] = serr
    return tab

def coverage_ratio(sample: np.ndarray, ref: np.ndarray) -> np.ndarray:
    # We compute the coverage ratio between sample and reference.
    # To estimate the standard error of the ratio, we use
    #
    #                     CV[X]^2 + CV[Y]^2 - CV[X]^2 CV[Y]^2
    #    CV[ X / Y ]^2 = -------------------------------------
    #                               1 - 2 * CV[Y]^2 
    #
    # where CV[X] = STDDEV[X] / MEAN[X] denotes the coefficient of variation
    # 
    # This identity holds exactly only for Gamma-distributed random
    # variables, but should hold approximately for any RV with short tails
    cols = [('q', np.float64), ('serr', np.float64)]
    tab = np.zeros(sample.shape[0], dtype=cols)
    tab["q"] = sample["avg"] / ref["avg"]
    s_cv2 = (sample["serr"] / sample["avg"]) ** 2
    r_cv2 = (ref["serr"] / ref["avg"]) ** 2
    q_cv2 = (s_cv2 + r_cv2 - s_cv2 * r_cv2) / (1 - 2 * r_cv2)
    tab["serr"] = np.sqrt(tab["q"] * q_cv2)
    return tab
    

def read_bam(bamfile):
    with pysam.AlignmentFile(bamfile, "rb") as bam:
        # Read sequence lengths from BAM header
        seqlengths = { seq: length for seq, length in zip(bam.header.references, bam.header.lengths) }
        # Create Coverage objects for all chromosomes
        cov = {seq: Coverage(length) for seq, length in seqlengths.items()}
        # Scan read pairs
        pairs = pyngstuff.read_pair_iterator(bam)
        for pair in pairs:
            if pair.forward.reference_start <= pair.reverse.reference_end:
                cov[pair.forward.reference_name].fragment_inclusive(
                  begin=pair.forward.reference_start, end=pair.reverse.reference_end - 1,
                  dir=1)
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

# Make output unbuffered
sys.stdout = io.TextIOWrapper(open(sys.stdout.fileno(), 'wb', 0), write_through=True)

# Command-line arguments
SAM_BAM_FILE = sys.argv[1]
SAM_MAPPABILITY_FILE = sys.argv[2]
REF_BAM_FILE = sys.argv[3]
REF_MAPPABILITY_FILE = sys.argv[4]
WINDOWSIZE = int(sys.argv[5])
ABUNDANCE_FILE = sys.argv[6]

# Read data
print(f"Reading sample BAM {SAM_BAM_FILE}")
SAM_COVERAGE, SAM_SEQLENGTHS = read_bam(SAM_BAM_FILE)

print(f"Reading reference BAM {REF_BAM_FILE}")
REF_COVERAGE, REF_SEQLENGTHS = read_bam(REF_BAM_FILE)

print(f"Reading mappabilities {SAM_MAPPABILITY_FILE}, {REF_MAPPABILITY_FILE}")
SAM_MAPPABILITY = read_mappability(SAM_MAPPABILITY_FILE, SAM_SEQLENGTHS)
REF_MAPPABILITY = read_mappability(REF_MAPPABILITY_FILE, REF_SEQLENGTHS)

# Compute
print(f"Finding mappable windows for window size {WINDOWSIZE}bp")

# NOTE: In the following, the data is kept in dictionaries which contain
# one numpy.ndarray per chromosome, indexed by the chromosome name. 

# Compute window midpoints
MIDPOINTS = { seq: np.arange(start=WINDOWSIZE/2, stop=length, step=WINDOWSIZE) for seq, length in SAM_SEQLENGTHS.items() }
# Test whether windows overlap some non-mappable regions
SAM_MAPPABLE = { seq: np.array([ not SAM_MAPPABILITY[seq].overlaps(i-WINDOWSIZE/2, i+WINDOWSIZE/2)
                                 for i in grid.tolist() ], dtype=np.bool_)
                 for seq, grid in MIDPOINTS.items() }
REF_MAPPABLE = { seq: np.array([ not REF_MAPPABILITY[seq].overlaps(i-WINDOWSIZE/2, i+WINDOWSIZE/2)
                                 for i in grid.tolist() ], dtype=np.bool_)
                 for seq, grid in MIDPOINTS.items() }

# Compute average coverage and std. dev. within each window
# window_average() returns a heterogeonous ndarray with named columns "avg", "sd" and "wlen"
print(f"Computing abundances and abundance ratios")
SAM = { seq: window_statistics(SAM_COVERAGE[seq], WINDOWSIZE)
        for seq, grid in MIDPOINTS.items() }
REF = { seq: window_statistics(REF_COVERAGE[seq], WINDOWSIZE)
        for seq, grid in MIDPOINTS.items() }

# Compute coverage ratio
RATIO = { seq: coverage_ratio(SAM[seq], REF[seq]) for seq in SAM.keys() }

# Output results
SLEN = max([len(seq) for seq, _ in MIDPOINTS.items() ])
COLS = [('seq', f"U{SLEN}"), ('position', np.int32), ('windowsize', np.int32), ('fragmentlength', np.float64), ('readcoverage', np.float64), 
        ("sample_avg", np.float64), ("sample_sd", np.float64), ("sample_avg_stderr", np.float64), ("sample_mappable", np.bool_), 
        ("reference_avg", np.float64), ("reference_sd", np.float64), ("reference_avg_stderr", np.float64), ("reference_mappable", np.bool_), 
        ("ratio", np.float64), ("ratio_stderr", np.float64), ("ratio_mappable", np.bool_) ]
with (open(ABUNDANCE_FILE, "w+") as f):
    f.write("\t".join([name for (name, type) in COLS]) + "\n")
    for seq in SAM.keys():
        tab = np.zeros(MIDPOINTS[seq].size, dtype=COLS)
        tab["seq"] = np.array([seq] * MIDPOINTS[seq].shape[0], dtype=tab["seq"].dtype)
        tab["position"] = MIDPOINTS[seq]
        tab["windowsize"] = SAM[seq]["wlen"]
        tab["fragmentlength"] = SAM[seq]["flen"]
        tab["readcoverage"] = SAM[seq]["fcov"]
        tab["sample_mappable"] = SAM_MAPPABLE[seq]
        tab["sample_avg"] = SAM[seq]["avg"]
        tab["sample_sd"] = SAM[seq]["sd"]
        tab["sample_avg_stderr"] = SAM[seq]["serr"]
        tab["reference_avg"] = REF[seq]["avg"]
        tab["reference_sd"] = REF[seq]["sd"]
        tab["reference_avg_stderr"] = REF[seq]["serr"]
        tab["reference_mappable"] = REF_MAPPABLE[seq]
        tab["ratio"] = RATIO[seq]["q"]
        tab["ratio_stderr"] = RATIO[seq]["serr"]
        tab["ratio_mappable"] = SAM_MAPPABLE[seq] & REF_MAPPABLE[seq]
        np.savetxt(f, tab, delimiter="\t", fmt=("%s",) + ("%d",)*2 + ("%.2f",)*2 + ("%.2f", "%.5f", "%.5f", "%d")*2 + ("%.5f", "%.5f", "%d"))
