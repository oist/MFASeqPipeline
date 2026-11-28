"""Tests that run mapping_statistics.py and abundance.py directly (without arubaito) on
small simulated BAM files, and compare their results with values computed here from
the simulated reads, following docs/stats_format.md and docs/abundance_format.md."""
import csv
import gzip
import math
import os
import statistics
import subprocess
import sys

import numpy as np
import pysam
import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAPPING_STATISTICS = os.path.join(REPO, "steps", "mapping_statistics.py")
ABUNDANCE = os.path.join(REPO, "steps", "abundance.py")

CHROMOSOMES = [("chrA", 5000), ("chrB", 3000)]
MAPQS = [0, 10, 20, 30, 60]

# Flags
PAIRED, PROPER, UNMAPPED, MATE_UNMAPPED, REVERSE, MATE_REVERSE = 0x1, 0x2, 0x4, 0x8, 0x10, 0x20
READ1, READ2, SECONDARY, SUPPLEMENTARY = 0x40, 0x80, 0x100, 0x800


def make_read(header, name, flag, chrom, start, length, mapq=60, mate_chrom=None, mate_start=0, tlen=0):
    r = pysam.AlignedSegment(header)
    r.query_name = name
    r.flag = flag
    r.reference_name = chrom
    r.reference_start = start
    r.mapping_quality = mapq
    r.query_sequence = "A" * length
    r.query_qualities = pysam.qualitystring_to_array("I" * length)
    if not flag & UNMAPPED:
        r.cigartuples = [(0, length)]
    if flag & PAIRED:
        r.next_reference_name = mate_chrom or chrom
        r.next_reference_start = mate_start
        r.template_length = tlen
    return r


def write_bam(path, chromosomes, make_reads):
    """Write a coordinate-sorted, indexed BAM with the reads returned by make_reads(header)."""
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"},
                                              "SQ": [{"SN": c, "LN": l} for c, l in chromosomes]})
    unsorted = str(path) + ".unsorted.bam"
    with pysam.AlignmentFile(unsorted, "wb", header=header) as bam:
        for r in make_reads(header):
            bam.write(r)
    pysam.sort("-o", str(path), unsorted)
    pysam.index(str(path))
    os.remove(unsorted)


def pair_reads(header, name, chrom, start, flen, rlen, mapq1, mapq2, read1_forward):
    """Reads of a properly paired fragment [start, start + flen), forward read first."""
    end = start + flen
    fwd_flag = PAIRED | PROPER | MATE_REVERSE | (READ1 if read1_forward else READ2)
    rev_flag = PAIRED | PROPER | REVERSE | (READ2 if read1_forward else READ1)
    return [make_read(header, name, fwd_flag, chrom, start, rlen, mapq1, chrom, end - rlen, flen),
            make_read(header, name, rev_flag, chrom, end - rlen, rlen, mapq2, chrom, start, -flen)]


# ---------------------------------------------------------------------------
# mapping_statistics.py

def test_mapping_statistics_paired(tmp_path):
    rng = np.random.RandomState(1)
    primary = []        # (is_mapped, is_proper, |tlen|, read length) of all primary alignments

    def reads(header):
        out = []
        for i in range(301):
            chrom, length = CHROMOSOMES[i % 2]
            flen, rlen = int(rng.randint(150, 400)), int(rng.choice([50, 75, 100]))
            start = int(rng.randint(0, length - flen))
            kind = i % 6
            if kind < 3:
                # Proper pair, one mate with a supplementary alignment
                pair = pair_reads(header, f"p{i}", chrom, start, flen, rlen, 60, 60, rng.rand() < 0.5)
                out += pair
                primary.extend([(True, True, flen, rlen), (True, True, flen, rlen)])
                if kind == 0:
                    out.append(make_read(header, f"p{i}", pair[0].flag | SUPPLEMENTARY, chrom, start, 20))
            elif kind == 3:
                # Pair that is not proper, plus a secondary alignment
                out.append(make_read(header, f"d{i}", PAIRED | READ1, chrom, start, rlen, 60, chrom, start + 1000, 0))
                out.append(make_read(header, f"d{i}", PAIRED | READ2, chrom, start + 1000, rlen, 60, chrom, start, 0))
                out.append(make_read(header, f"d{i}", PAIRED | READ1 | SECONDARY, chrom, start + 50, rlen))
                primary.extend([(True, False, 0, rlen), (True, False, 0, rlen)])
            elif kind == 4:
                # One mate unmapped, placed at the position of the other
                out.append(make_read(header, f"u{i}", PAIRED | READ1 | MATE_UNMAPPED, chrom, start, rlen, 60,
                                     chrom, start, 0))
                out.append(make_read(header, f"u{i}", PAIRED | READ2 | UNMAPPED, chrom, start, rlen, 0,
                                     chrom, start, 0))
                primary.extend([(True, False, 0, rlen), (False, False, 0, rlen)])
            else:
                # Proper pair with low mapping quality, which counts like any other
                out += pair_reads(header, f"q{i}", chrom, start, flen, rlen, 0, 3, rng.rand() < 0.5)
                primary.extend([(True, True, flen, rlen), (True, True, flen, rlen)])
        return out

    write_bam(tmp_path / "s.bam", CHROMOSOMES, reads)
    subprocess.run([sys.executable, MAPPING_STATISTICS, tmp_path / "s.bam", tmp_path / "s.stats"], check=True)

    with open(tmp_path / "s.stats") as f:
        stats = dict(line.rstrip("\n").split("\t") for line in f)
    assert list(stats) == ["total_reads", "mapped_reads", "properly_paired_reads", "insert_size", "read_length"]
    assert int(stats["total_reads"]) == len(primary)
    assert int(stats["mapped_reads"]) == sum(1 for m, _, _, _ in primary if m)
    assert int(stats["properly_paired_reads"]) == sum(1 for _, p, _, _ in primary if p)
    assert int(stats["insert_size"]) == int(statistics.median([t for _, p, t, _ in primary if p]))
    assert int(stats["read_length"]) == int(statistics.median([l for _, _, _, l in primary]))


def test_mapping_statistics_single(tmp_path):
    rng = np.random.RandomState(2)
    lengths = []

    def reads(header):
        out = []
        for i in range(200):
            chrom, length = CHROMOSOMES[i % 2]
            rlen = int(rng.randint(30, 150))
            start = int(rng.randint(0, length - rlen))
            out.append(make_read(header, f"r{i}", REVERSE if i % 3 == 0 else 0, chrom, start, rlen))
            lengths.append(rlen)
        out.append(make_read(header, "r0", SECONDARY, "chrA", 10, 40))
        out.append(make_read(header, "unmapped", UNMAPPED, "chrA", 10, 60))
        lengths.append(60)
        return out

    write_bam(tmp_path / "s.bam", CHROMOSOMES, reads)
    subprocess.run([sys.executable, MAPPING_STATISTICS, tmp_path / "s.bam", tmp_path / "s.stats"], check=True)

    with open(tmp_path / "s.stats") as f:
        stats = dict(line.rstrip("\n").split("\t") for line in f)
    assert stats == {"total_reads": "201", "mapped_reads": "200", "properly_paired_reads": "0",
                     "insert_size": "nan", "read_length": str(int(statistics.median(lengths)))}


# ---------------------------------------------------------------------------
# abundance.py

def simulate_fragments(rng, n, chromosomes, paired):
    """Random fragments (chrom, start, flen, rlen, mapq1, mapq2, read1_forward) along the chromosomes."""
    frags = []
    for i in range(n):
        chrom, length = chromosomes[i % len(chromosomes)]
        rlen = int(rng.choice([50, 100]))
        flen = int(rng.randint(rlen, 400)) if paired else rlen
        start = int(rng.randint(0, length - flen + 1))
        mapq1 = int(rng.choice(MAPQS))
        # Single reads only have one mapping quality
        mapq2 = int(rng.choice(MAPQS)) if paired else mapq1
        frags.append((chrom, start, flen, rlen, mapq1, mapq2, bool(rng.rand() < 0.5)))
    return frags


def write_fragments_bam(path, chromosomes, frags, paired, rng):
    """Write the fragments as read pairs or single reads, plus alignments abundance.py must ignore."""
    def reads(header):
        out = []
        for i, (chrom, start, flen, rlen, mapq1, mapq2, read1_forward) in enumerate(frags):
            if paired:
                out += pair_reads(header, f"f{i}", chrom, start, flen, rlen, mapq1, mapq2, read1_forward)
            else:
                out.append(make_read(header, f"f{i}", 0 if read1_forward else REVERSE, chrom, start, rlen, mapq1))
        for i in range(50):
            chrom, length = chromosomes[i % len(chromosomes)]
            start = int(rng.randint(0, length - 500))
            # Secondary and supplementary alignments, unmapped reads
            out.append(make_read(header, f"f{i}", SECONDARY | (PAIRED | PROPER | READ1 if paired else 0),
                                 chrom, start, 50, 60, chrom, start + 200, 250))
            out.append(make_read(header, f"f{i}", SUPPLEMENTARY | (PAIRED | PROPER | READ1 if paired else 0),
                                 chrom, start, 50, 60, chrom, start + 200, 250))
            out.append(make_read(header, f"x{i}", UNMAPPED, chrom, start, 50, 0))
            if paired:
                # Pairs that are not proper, mapped on different chromosomes, or with one mate unmapped
                out.append(make_read(header, f"d{i}", PAIRED | READ1, chrom, start, 50, 60, chrom, start + 300, 0))
                out.append(make_read(header, f"d{i}", PAIRED | READ2 | REVERSE, chrom, start + 300, 50, 60,
                                     chrom, start, 0))
                other = chromosomes[(i + 1) % len(chromosomes)][0]
                out.append(make_read(header, f"c{i}", PAIRED | PROPER | READ1, chrom, start, 50, 60, other, 100, 0))
                out.append(make_read(header, f"m{i}", PAIRED | READ1 | MATE_UNMAPPED, chrom, start, 50, 60,
                                     chrom, start, 0))
        return out
    write_bam(path, chromosomes, reads)


def write_mappability(path, chromosomes, regions):
    with open(path, "w") as f:
        print("# Command used to create this file:", file=f)
        print("#   simulated", file=f)
        print("# Chromosomes:", file=f)
        for chrom, length in chromosomes:
            print(f"#   {chrom}: {length}", file=f)
        for chrom, start, end in regions:
            print(f"{chrom}\t{start}\t{end}", file=f)


def expected_windows(chromosomes, frags, mapq, windowsize, regions):
    """Per-window statistics as described in docs/abundance_format.md, computed from the fragments."""
    result = {}
    for chrom, length in chromosomes:
        counted = [(s, s + f - 1) for (c, s, f, r, q1, q2, _) in frags if c == chrom and q1 >= mapq and q2 >= mapq]
        cov = np.zeros(length)
        sp = np.zeros(length, dtype=int)
        ep = np.zeros(length, dtype=int)
        for first, last in counted:
            cov[first:last + 1] += 1
            sp[first] += 1
            ep[last] += 1
        flen = round(sum(l - f + 1 for f, l in counted) / len(counted))
        rows = []
        for first in range(0, length, windowsize):
            last = min(first + windowsize, length) - 1
            n = last - first + 1
            # Mean number of bases of a fragment of length flen in the window, over all
            # placements of the fragment that overlap the window
            overlaps = [min(i + flen, n) - max(i, 0) for i in range(-flen + 1, n)]
            fcov = sum(overlaps) / len(overlaps)
            avg = cov[first:last + 1].mean()
            rows.append({
                "seq": chrom, "first": first, "midpoint": (first + last) / 2, "last": last, "windowsize": n,
                "fragmentlength": flen, "fragmentcoverage": fcov,
                "coverage_avg": avg, "coverage_sd": cov[first:last + 1].std(),
                "coverage_avg_stderr": math.sqrt(avg) * math.sqrt(fcov / n),
                "startingpoints": int(sp[first:last + 1].sum()), "endpoints": int(ep[first:last + 1].sum()),
                "mappable": int(not any(c == chrom and s <= last and e > first for c, s, e in regions)),
            })
        result[chrom] = rows
    return result


def cv2(stderr, mean):
    """Squared coefficient of variation stderr / mean (infinite for mean zero, nan if both are zero)."""
    if mean == 0:
        return math.nan if stderr == 0 else math.inf
    return (stderr / mean) ** 2


def expected_ratio(x, y, x_cv2, y_cv2):
    """Ratio x / y and its standard error, as in docs/abundance_format.md (nan if undefined)."""
    if y == 0:
        return math.nan, math.nan
    q = x / y
    if not 1 - 2 * y_cv2 > 0:
        return q, math.nan
    q_cv2 = (x_cv2 + y_cv2 - x_cv2 * y_cv2) / (1 - 2 * y_cv2)
    if not math.isfinite(q_cv2):
        return q, math.nan
    return q, q * math.sqrt(q_cv2)


def read_abundance(path):
    with gzip.open(path, "rt") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def assert_value(row, column, expected, decimals):
    actual = float(row[column])
    if isinstance(expected, float) and math.isnan(expected):
        assert math.isnan(actual), f"{column} at {row['seq']}:{row['first']} is {actual}, expected nan"
    else:
        assert abs(actual - expected) <= 0.5 * 10 ** -decimals + 1e-9, \
            f"{column} at {row['seq']}:{row['first']} is {actual}, expected {expected}"


def assert_sample_columns(row, exp, prefix):
    assert_value(row, f"{prefix}_coverage_avg", exp["coverage_avg"], 2)
    assert_value(row, f"{prefix}_coverage_sd", exp["coverage_sd"], 5)
    assert_value(row, f"{prefix}_coverage_avg_stderr", exp["coverage_avg_stderr"], 5)
    assert int(row[f"{prefix}_startingpoints"]) == exp["startingpoints"]
    assert int(row[f"{prefix}_endpoints"]) == exp["endpoints"]
    assert int(row[f"{prefix}_mappable"]) == exp["mappable"]


def check_abundance(path, chromosomes, sample, reference):
    rows = read_abundance(path)
    assert len(rows) == sum(len(sample[c]) for c, _ in chromosomes)
    expected_rows = [(e, reference[c][i] if reference else None)
                     for c, _ in chromosomes for i, e in enumerate(sample[c])]
    for row, (s, r) in zip(rows, expected_rows):
        assert row["seq"] == s["seq"]
        assert int(row["first"]) == s["first"]
        assert int(row["last"]) == s["last"]
        assert int(row["windowsize"]) == s["windowsize"]
        assert_value(row, "midpoint", s["midpoint"], 1)
        assert_value(row, "fragmentlength", s["fragmentlength"], 2)
        assert_value(row, "fragmentcoverage", s["fragmentcoverage"], 2)
        assert_sample_columns(row, s, "sample")
        if r is None:
            for column in ["reference_coverage_avg", "reference_coverage_sd", "reference_coverage_avg_stderr",
                           "coverage_ratio", "coverage_ratio_stderr", "startingpoint_ratio",
                           "startingpoint_ratio_stderr", "endpoint_ratio", "endpoint_ratio_stderr"]:
                assert_value(row, column, math.nan, 5)
            assert int(row["reference_startingpoints"]) == -1
            assert int(row["reference_endpoints"]) == -1
            assert int(row["reference_mappable"]) == 0
            assert int(row["ratio_mappable"]) == 0
            continue
        assert_sample_columns(row, r, "reference")
        q, q_se = expected_ratio(s["coverage_avg"], r["coverage_avg"],
                                 cv2(s["coverage_avg_stderr"], s["coverage_avg"]),
                                 cv2(r["coverage_avg_stderr"], r["coverage_avg"]))
        assert_value(row, "coverage_ratio", q, 5)
        assert_value(row, "coverage_ratio_stderr", q_se, 5)
        # For counts, CV^2 = 1 / count
        q, q_se = expected_ratio(s["startingpoints"], r["startingpoints"],
                                 cv2(1, s["startingpoints"]) ** 0.5, cv2(1, r["startingpoints"]) ** 0.5)
        assert_value(row, "startingpoint_ratio", q, 5)
        assert_value(row, "startingpoint_ratio_stderr", q_se, 5)
        q, q_se = expected_ratio(s["endpoints"], r["endpoints"],
                                 cv2(1, s["endpoints"]) ** 0.5, cv2(1, r["endpoints"]) ** 0.5)
        assert_value(row, "endpoint_ratio", q, 5)
        assert_value(row, "endpoint_ratio_stderr", q_se, 5)
        assert int(row["ratio_mappable"]) == (s["mappable"] and r["mappable"])


def test_abundance_paired_with_reference(tmp_path):
    rng = np.random.RandomState(3)
    # The sample has twice the coverage of the reference near the start of each chromosome
    sample = simulate_fragments(rng, 3000, CHROMOSOMES, paired=True)
    sample += simulate_fragments(rng, 1000, [("chrA", 2000), ("chrB", 1200)], paired=True)
    reference = simulate_fragments(rng, 3000, CHROMOSOMES, paired=True)
    sample_regions = [("chrA", 1000, 1200), ("chrB", 2990, 3000)]
    reference_regions = [("chrA", 4100, 4101), ("chrA", 1150, 1300)]
    write_fragments_bam(tmp_path / "s.bam", CHROMOSOMES, sample, True, rng)
    write_fragments_bam(tmp_path / "r.bam", CHROMOSOMES, reference, True, rng)
    write_mappability(tmp_path / "s.bed", CHROMOSOMES, sample_regions)
    write_mappability(tmp_path / "r.bed", CHROMOSOMES, reference_regions)

    subprocess.run([sys.executable, ABUNDANCE, "--sample-bam", tmp_path / "s.bam",
                    "--sample-mappability", tmp_path / "s.bed", "--reference-bam", tmp_path / "r.bam",
                    "--reference-mappability", tmp_path / "r.bed", "--paired", "--mapping-quality", "20",
                    "--windowsize", "700", "--output", tmp_path / "s.abd.gz"], check=True)

    check_abundance(tmp_path / "s.abd.gz", CHROMOSOMES,
                    expected_windows(CHROMOSOMES, sample, 20, 700, sample_regions),
                    expected_windows(CHROMOSOMES, reference, 20, 700, reference_regions))


def test_abundance_paired_without_reference(tmp_path):
    rng = np.random.RandomState(4)
    sample = simulate_fragments(rng, 3000, CHROMOSOMES, paired=True)
    regions = [("chrB", 0, 10), ("chrB", 1500, 1501)]
    write_fragments_bam(tmp_path / "s.bam", CHROMOSOMES, sample, True, rng)
    write_mappability(tmp_path / "s.bed", CHROMOSOMES, regions)

    subprocess.run([sys.executable, ABUNDANCE, "--sample-bam", tmp_path / "s.bam",
                    "--sample-mappability", tmp_path / "s.bed", "--paired", "--mapping-quality", "30",
                    "--windowsize", "1k", "--output", tmp_path / "s.abd.gz"], check=True)

    check_abundance(tmp_path / "s.abd.gz", CHROMOSOMES, expected_windows(CHROMOSOMES, sample, 30, 1000, regions), None)


def test_abundance_single_with_reference(tmp_path):
    rng = np.random.RandomState(5)
    sample = simulate_fragments(rng, 4000, CHROMOSOMES, paired=False)
    reference = simulate_fragments(rng, 3000, CHROMOSOMES, paired=False)
    regions = [("chrA", 0, 5000)]
    write_fragments_bam(tmp_path / "s.bam", CHROMOSOMES, sample, False, rng)
    write_fragments_bam(tmp_path / "r.bam", CHROMOSOMES, reference, False, rng)
    write_mappability(tmp_path / "s.bed", CHROMOSOMES, [])
    write_mappability(tmp_path / "r.bed", CHROMOSOMES, regions)

    subprocess.run([sys.executable, ABUNDANCE, "--sample-bam", tmp_path / "s.bam",
                    "--sample-mappability", tmp_path / "s.bed", "--reference-bam", tmp_path / "r.bam",
                    "--reference-mappability", tmp_path / "r.bed", "--mapping-quality", "10",
                    "--windowsize", "600", "--output", tmp_path / "s.abd.gz"], check=True)

    check_abundance(tmp_path / "s.abd.gz", CHROMOSOMES,
                    expected_windows(CHROMOSOMES, sample, 10, 600, []),
                    expected_windows(CHROMOSOMES, reference, 10, 600, regions))


def test_abundance_uncompressed_output(tmp_path):
    rng = np.random.RandomState(6)
    sample = simulate_fragments(rng, 500, CHROMOSOMES, paired=True)
    write_fragments_bam(tmp_path / "s.bam", CHROMOSOMES, sample, True, rng)
    write_mappability(tmp_path / "s.bed", CHROMOSOMES, [])

    subprocess.run([sys.executable, ABUNDANCE, "--sample-bam", tmp_path / "s.bam",
                    "--sample-mappability", tmp_path / "s.bed", "--paired", "--mapping-quality", "0",
                    "--windowsize", "5k", "--output", tmp_path / "s.abd"], check=True)

    with open(tmp_path / "s.abd") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    assert [(r["seq"], r["first"], r["last"]) for r in rows] == [("chrA", "0", "4999"), ("chrB", "0", "2999")]


def test_abundance_mismatched_chromosomes(tmp_path):
    rng = np.random.RandomState(7)
    other = [("chrA", 5000), ("chrB", 3100), ("chrC", 100)]
    write_fragments_bam(tmp_path / "s.bam", CHROMOSOMES, simulate_fragments(rng, 100, CHROMOSOMES, True), True, rng)
    write_bam(tmp_path / "r.bam", other, lambda header: pair_reads(header, "p", "chrA", 100, 300, 50, 60, 60, True))
    write_mappability(tmp_path / "s.bed", CHROMOSOMES, [])
    write_mappability(tmp_path / "r.bed", other, [])

    r = subprocess.run([sys.executable, ABUNDANCE, "--sample-bam", tmp_path / "s.bam",
                        "--sample-mappability", tmp_path / "s.bed", "--reference-bam", tmp_path / "r.bam",
                        "--reference-mappability", tmp_path / "r.bed", "--paired", "--mapping-quality", "20",
                        "--windowsize", "1k", "--output", tmp_path / "s.abd.gz"],
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    assert r.returncode == 1
    assert "chrC (100 bp) is missing from the sample BAM" in r.stderr
    assert "chrB has length 3000 bp in the sample BAM but 3100 bp in the reference BAM" in r.stderr
    assert not os.path.exists(tmp_path / "s.abd.gz")


@pytest.mark.parametrize("args, message", [
    (["--reference-bam", "r.bam"], "--reference-bam and --reference-mappability must be specified together"),
    (["--windowsize", "10x"], "invalid window size '10x'"),
])
def test_abundance_invalid_arguments(tmp_path, args, message):
    base = {"--sample-bam": "s.bam", "--sample-mappability": "s.bed", "--mapping-quality": "20",
            "--windowsize": "1k", "--output": "s.abd.gz"}
    for i in range(0, len(args), 2):
        base[args[i]] = args[i + 1]
    r = subprocess.run([sys.executable, ABUNDANCE] + [x for kv in base.items() for x in kv], cwd=tmp_path,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    assert r.returncode == 2
    assert message in r.stderr
