import sys
import math

import numpy
import pyngstuff
import pysam
from intervaltree import IntervalTree, Interval

class Abundance:
    def __init__(self, length, resolution):
        self.length = length
        self.resolution = resolution
        # The positions at which we compute abundances
        self.positions = [ (i-1)*resolution for i in range(1, math.floor(length / resolution)) ]
        
        # The read pairs that begin ]
        self.begin = []
        self.abundance = [math.nan] * len(self.postions)
        self.current = 0
        self.begin = []
        self.end = []

    def increment_range(self, begin, end):
        assert (len(self.begin) == 0) or self.begin[0][0] <= begin
        assert begin <= end
        self.begin.append(begin)
        self.end.append(end)

    def summarize(self, upto):
        while (len(self.begin) > 0) and self.begin[0] <= upto:
            self.begin.


BAMFILE = sys.argv[1]
MAPPABILITYFILE = sys.argv[2]
ABUNDANCEFILE = sys.argv[3]

# Open BAM file
bam = pysam.AlignmentFile(BAMFILE, "rb")

# Read sequence lengths from header
SEQLENGTHS = [ {seq: length} for seq, length in zip(bam.header.references, bam.header.lengths)]

# Read mappability data
MAPPABILITY = [ {seq: IntervalTree([Interval(0, length)]) } for seq, length in SEQLENGTHS.items() ]
with open(MAPPABILITYFILE, "r") as f:
    for line in f.readlines():
        (seq, start, end) = line.split("\t")
        MAPPABILITY[seq].chop(start, end)

# Scan read pairs
read_pairs = pyngstuff.read_pair_ordered_iterator(pyngstuff.read_pair_iterator(bam))
startingpoints =
for pair in read_pairs:
