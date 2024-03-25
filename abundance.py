import sys
import numpy
import pyngstuff
import pysam
from intervaltree import IntervalTree, Interval

class Coverage:
    def __init__(self, length : int):
        self.length = length
        self.deltas = numpy.zeros([length], dtype=numpy.int32)

    def range_inclusive(self, begin : int, end : int, dir : int):
        assert (begin >= 0) and (begin < self.length)
        assert (end >= 0) and (end < self.length)
        assert dir in [-1, 1]
        if dir == -1:
            self.range_inclusive(end, begin, +1)
        if begin <= end:
            # [begin..end]
            self.deltas[begin] += 1
            if end < self.length - 1:
                self.deltas[end+1] -= 1
        elif (begin > end) and (dir == 1):
            # [0..end] u [begin...last]
            self.deltas[0] += 1
            self.deltas[end+1] -= -1
            self.deltas[begin] += 1

    def coverage(self, positions: numpy.ndarray) -> numpy.ndarray:
        cov = numpy.cumsum(self.deltas)
        return cov[positions]


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
                cov[pair.forward.reference_name].range_inclusive(begin=pair.forward.reference_start,
                                                                 end=pair.reverse.reference_end - 1,
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

SAMPLE_BAM_FILE = sys.argv[1]
SAMPLE_MAPPABILITY_FILE = sys.argv[2]
REFERENCE_BAM_FILE = sys.argv[3]
REFERENCE_MAPPABILITY_FILE = sys.argv[4]
RESOLUTION = int(sys.argv[5])
PADDING = int(sys.argv[6])
ABUNDANCE_FILE = sys.argv[7]

print(f"Reading sample BAM {SAMPLE_BAM_FILE}")
SAMPLE_COVERAGE, SAMPLE_SEQLENGTHS = read_bam(SAMPLE_BAM_FILE)

print(f"Reading reference BAM {REFERENCE_BAM_FILE}")
REFERENCE_COVERAGE, REFERENCE_SEQLENGTHS = read_bam(REFERENCE_BAM_FILE)

print(f"Reading mappabilities {SAMPLE_MAPPABILITY_FILE}, {REFERENCE_MAPPABILITY_FILE}")
SAMPLE_MAPPABILITY = read_mappability(SAMPLE_MAPPABILITY_FILE, SAMPLE_SEQLENGTHS)
REFERENCE_MAPPABILITY = read_mappability(REFERENCE_MAPPABILITY_FILE, REFERENCE_SEQLENGTHS)

print(f"Finding mappable grid points for grid size {RESOLUTION}bp")
GRID = { seq: numpy.arange(start=0, stop=length, step=RESOLUTION) for seq, length in SAMPLE_SEQLENGTHS.items() }
SAMPLE_MAPPABLE = { seq: numpy.array([ not SAMPLE_MAPPABILITY[seq].overlaps(i-PADDING, i+PADDING+1)
                                       for i in grid.tolist() ], dtype=numpy.bool_)
                    for seq, grid in GRID.items() }
REFERENCE_MAPPABLE = { seq: numpy.array([ not REFERENCE_MAPPABILITY[seq].overlaps(i-PADDING, i+PADDING+1)
                                       for i in grid.tolist() ], dtype=numpy.bool_)
                    for seq, grid in GRID.items() }

print(f"Computing abundances and abundance ratios")
SAMPLE = { seq: numpy.where(SAMPLE_MAPPABLE[seq], SAMPLE_COVERAGE[seq].coverage(grid), numpy.nan)
           for seq, grid in GRID.items() }
REF = { seq: numpy.where(REFERENCE_MAPPABLE[seq], REFERENCE_COVERAGE[seq].coverage(grid), numpy.nan)
        for seq, grid in GRID.items() }
RATIO = { seq: SAMPLE[seq] / REF[seq] for seq in SAMPLE.keys() }
STDDEV = { seq: ratios * numpy.sqrt(1.0 / SAMPLE[seq] + 1.0 / REF[seq]) for seq, ratios in RATIO.items() }

SLEN = max([len(seq) for seq, _ in GRID.items() ])
COLS = [('seq', f"U{SLEN}"), ('pos', numpy.int32), ("sample", numpy.float64), ("reference", numpy.float64),
        ("ratio", numpy.float64), ("stddev", numpy.float64)]
with (open(ABUNDANCE_FILE, "w+") as f):
    f.write("\t".join([name for (name, type) in COLS]) + "\n")
    for seq in SAMPLE.keys():
        tab = numpy.zeros(GRID[seq].size, dtype=COLS)
        tab["seq"] = numpy.array([seq] * GRID[seq].shape[0], dtype=tab["seq"].dtype)
        tab["pos"] = GRID[seq]
        tab["sample"] = SAMPLE[seq]
        tab["reference"] = REF[seq]
        tab["ratio"] = RATIO[seq]
        tab["stddev"] = STDDEV[seq]
        numpy.savetxt(f, tab, delimiter="\t", fmt=("%s", "%d", "%.0f", "%.0f", "%.5f", "%.5f"))