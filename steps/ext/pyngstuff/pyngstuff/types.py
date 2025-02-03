from collections import namedtuple
from enum import Enum

class Strand(Enum):
    PLUS = 1
    MINUS = 2

    def __str__(self):
        if self.value == 1:
            return "+"
        else:
            return "-"

class UmiFuzzy(namedtuple('UmiFuzzy', ['pos1', 'pos2', 'tag'])):
    def __new__(cls, pos1=-1, pos2=-1, tag=""):
        return super().__new__(cls,
                               pos1 if pos1 is not None else -1,
                               pos2 if pos2 is not None else -1,
                               tag if tag is not None else "")

class UmiFixed(namedtuple('UmiFixed', ['sample', 'seq', 'strand'])):
    def __new__(cls, sample=None, seq=None, strand=None):
        return super().__new__(cls, sample, seq, strand)

class Umi(namedtuple('Umi', ['sample', 'seq', 'strand', 'pos1', 'pos2', 'tag'])):
    """ The full UMI """
    def __new__(cls, sample=None, seq=None, strand=None, pos1=-1, pos2=-1, tag=""):
        return super().__new__(cls, sample, seq, strand,
                               pos1 if pos1 is not None else -1,
                               pos2 if pos2 is not None else -1,
                               tag if tag is not None else "")

    @staticmethod
    def combine(fixed, fuzzy):
        return Umi(fixed.sample, fixed.seq, fixed.strand, fuzzy.pos1, fuzzy.pos2, fuzzy.tag)

    @property
    def fuzzy(self):
        return UmiFuzzy(self.pos1, self.pos2, self.tag)

    @property
    def fixed(self):
        return UmiFixed(self.sample, self.seq, self.strand)

    def __str__(self):
        return ''.join([
            "%s#" % self.sample if self.sample is not None else "",
            "%s:" % self.seq if self.seq is not None else "",
            "%d" % self.pos1 if self.pos1 >= 0 else "",
            "-%d" % self.pos2 if self.pos2 >= 0 else "",
            "(%s)" % self.strand if self.strand is not None else "",
            "|%s" % self.tag if self.tag != "" else ""
        ])

    def __repr__(self):
        return str(self)
