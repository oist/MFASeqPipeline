import abc
import regex
import argparse
import sys
from collections import namedtuple, defaultdict

class UmiTagSourceAutodetectionFailure(Exception):
    def __init__(self):
        super().__init__("autodetection of the UMI tag source failed, please specify manually")

class UmiTagSource:
    @staticmethod
    def parse(spec):
        if spec == "auto":
            return UmiTagSourceAutodetect()
        elif spec[:8] == "readname":
            return UmiTagSourceReadname.parse(spec)
        else:
            return UmiTagSourceBAMTags.parse(spec)

    @staticmethod
    def autodetect(read, cellbarcode=False):
        source = UmiTagSourceBAMTags.autodetect(read, cellbarcode)
        if source is not None:
            return source
        source = UmiTagSourceReadname.autodetect(read, cellbarcode)
        if source is not None:
            return source
        raise UmiTagSourceAutodetectionFailure()

    @property
    @abc.abstractmethod
    def cellbarcode(self):
        pass


class UmiTagSourceAutodetect(UmiTagSource):
    def __init__(self):
        self._cellbarcode = None
        self.source = None

    @property
    def cellbarcode(self):
        return self._cellbarcode

    @cellbarcode.setter
    def cellbarcode(self, value):
        if self.source is not None:
            return RuntimeError("cannot change cellbarcode setting after auto-detection has taken place")
        self._cellbarcode = value

    def __call__(self, read):
        if self.source is None:
            print("*** Auto-detecting the UMI tag source... ", file=sys.stdout, end='', flush=True)
            if self.cellbarcode is None:
                try:
                    self.source = UmiTagSource.autodetect(read, True)
                except UmiTagSourceAutodetectionFailure:
                    self.cellbarcode = False
            if self.source is None:
                self.source = UmiTagSource.autodetect(read, self.cellbarcode)
            print("wil use %s" % str(self.source), file=sys.stdout, flush=True)
        return self.source(read)


class UmiTagSourceBAMTags(UmiTagSource, namedtuple('UmiTagSourceBAMTags', ['bamtag_umitag', 'bamtag_cellbc'])):
    @staticmethod
    def parse(spec):
        # Parse tag source specified, which has the form <var>=<tag>,... where
        # var is either C or U, and tag the BAM tag to read the UMI tag (U) or
        # cell barcode (C) from.
        m = regex.match(r'(?:(?P<var>[CU])=(?P<tag>[A-Za-z][A-Za-z0-9]),)*(?P<var>[CU])=(?P<tag>[A-Za-z][A-Za-z0-9])(?:,|$)+', spec)
        if m is None:
            raise argparse.ArgumentTypeError("invalid UMI tag source %s, expected a specification like U=RX,C=CB" % spec)
        if len(m.captures("var")) != len(set(m.captures("var"))):
            raise argparse.ArgumentTypeError("invalid UMI tag source %s, each of U and C may appear only once left of '='" % spec)

        # Create UmiTagSourceBAMTags instance
        bamtags = defaultdict(lambda: None)
        for var, tag in zip(m.captures("var"), m.captures("tag")):
            bamtags[var] = tag
        return UmiTagSourceBAMTags(bamtags['U'], bamtags['C'])

    @staticmethod
    def autodetect(read, cellbarcode=False):
        # Try known BAM tags that may contain UMI tags
        bamtag_umitag = None
        for bamtag in ["RX", "OX"]:
            try:
                v = read.get_tag(bamtag)
                if isinstance(v, str) and regex.match("^[ACGTNacgtn-]*$", v):
                    bamtag_umitag = bamtag
                    break
            except KeyError:
                pass
        if not cellbarcode:
            return UmiTagSourceBAMTags(bamtag_umitag, None) \
                if bamtag_umitag is not None else None

        # Try known BAM tags that may contain cell barcodes
        bamtag_cellbc = None
        for bamtag in ["CB", "CR"]:
            try:
                v = read.get_tag(bamtag)
                if isinstance(v, str) and regex.match("^[ACGTNacgtn-]*$", v):
                    bamtag_cellbc= bamtag
                    break
            except KeyError:
                pass

        # If we found all the necessary tags, we succeeded
        return UmiTagSourceBAMTags(bamtag_umitag, bamtag_cellbc) \
            if bamtag_umitag is not None and bamtag_cellbc is not None else None


    def __call__(self, read):
        return (read.get_tag(self.bamtag_umitag) if self.bamtag_umitag is not None else None,
                read.get_tag(self.bamtag_cellbc) if self.bamtag_cellbc is not None else None)

    def __str__(self):
        return ",".join("%s=%s" % (var, tag) for var, tag in (("U", self.bamtag_umitag), ("C", self.bamtag_cellbc))
                        if tag is not None)

    @property
    def cellbarcode(self):
        return self.bamtag_cellbc is not None


class UmiTagSourceReadname(UmiTagSource, namedtuple('UmiTagSourceReadname',
                                                    ['pattern', 'regexp', 'has_umitag', 'has_cellbc'])):
    @staticmethod
    def parse(pattern):
        # Parse readname pattern, which has the form readname<sep><var>...,
        # where var is either C (cell barcode) or U (UMI tag), and sep a single character.
        m = regex.match(r'^readname(?:(.)([CU]))+$', pattern)
        if m is None:
            raise argparse.ArgumentTypeError("invalid readname pattern %s, expected a pattern like readname_C_U" % pattern)

        regexp = UmiTagSourceReadname.make_regex((sep, var) for sep, var in zip(m.captures(1), m.captures(2)))
        return UmiTagSourceReadname(pattern, regexp, "U" in m.captures(2), "C" in m.captures(2))

    @staticmethod
    def make_regex(spec):
        # Translate pattern into a regular expression. Since the same var may
        # appear multiple times, we use the regex module, not the built-in re module,
        # since it allows us to access *all* matches of a group
        return regex.compile(r'^[^\s]+?%s(?:\s+.*)?$' % "".join(r'%s(?P<%s>[ACGTNacgtn-]+)' % (regex.escape(sep), var)
                                                                for sep, var in spec))

    @staticmethod
    def autodetect(read, cellbarcodes=False):
        # Assume the readname has the form "...<sep><cell barcode><sep><umi tag> ..."
        if not cellbarcodes:
            m = regex.match(r'^[^\s]+?(?P<U>(?:[^\w]|_))([ACGTNacgtn]+-)*[ACGTNacgtn]+(\s+.*)?$',
                            read.query_name, regex.ASCII)
            if m is None:
                return None
            pattern = "readname%sU" % m.group("U")
            regexp = UmiTagSourceReadname.make_regex((sep, var) for sep, var in ((m.group("U"), "U"),))
            return UmiTagSourceReadname(pattern, regexp, True, False)
        else:
            # Note that non-greed match operator "+?", which ensured that we try to
            # match as much as possible of the readname as a cell barcode or UMI tag.
            m = regex.match(r'^[^\s]+?(?P<C>(?:[^\w]|_))([ACGTNacgtn]+-)*[ACGTNacgtn]+(?P<U>(?:[^\w-]|_))([ACGTNacgtn]+-)*[ACGTNacgtn]+(\s+.*)?$',
                            read.query_name, regex.ASCII)
            if m is None:
                return None
            pattern = "readname%sC%sU" % (m.group("C"), m.group("U"))
            regexp = UmiTagSourceReadname.make_regex((sep, var) for sep, var in ((m.group("C"), "C"), (m.group("U"), "U")))
            return UmiTagSourceReadname(pattern, regexp, True, True)

    def __call__(self, read):
        m = self.regexp.match(read.query_name)
        if m is None:
            raise RuntimeError("pattern %s does not match read name %s" % (self.pattern, read.query_name))
        return ("-".join(m.captures("U")) if self.has_umitag else None,
                "-".join(m.captures("C")) if self.has_cellbc else None)

    def __str__(self):
        return self.pattern

    @property
    def cellbarcode(self):
        return self.has_cellbc

