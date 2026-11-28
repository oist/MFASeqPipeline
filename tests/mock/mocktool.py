#!/usr/bin/env python3
"""Mock for the tools called by the pipeline steps.

Each mock is installed under the name of the tool it replaces (a symlink to this
file in tests/mock/bin, or a file of the same name in the mock pipeline directory
tests/mock/pipeline), and figures out which tool to imitate from that name.

A mock does not do any real work. Instead, it reads the seed from the environment
variable MOCK_SEED, hashes all input files with an HMAC keyed by the seed, and
writes each output file as a "record" listing the tool, its arguments, the seed
and the input hashes. Records are text files with one tab-separated entry per line:

    tool    <name>
    arg     <argument>                  (one line per argument, in order)
    seed    <seed>
    input   <path>  <hash>              (one line per input file)
    stdin   <hash>                      (if the tool read its standard input)
    | <line>                            (the standard input, verbatim)
    <key>   <value>                     (tool-specific extra lines)

Input files are recognized by being existing regular files or named pipes given
as arguments (or as the value of --option=value). Compressed (gzip, bzip2) files
are hashed after decompression. Output files are written compressed if their name
ends in .gz or .bz2.
"""
import hashlib
import hmac
import os
import re
import stat
import sys


def fail(msg):
    print(f"MOCK ERROR ({TOOL}): {msg}", file=sys.stderr)
    sys.exit(2)


# gzip, bz2, glob and zipfile are imported only where needed, which saves noticeable time
# on some systems. Mocks are started many times by each test.

def decompress(data):
    if data.startswith(b"\x1f\x8b"):
        import gzip
        return gzip.decompress(data)
    if data.startswith(b"BZh"):
        import bz2
        return bz2.decompress(data)
    return data


def seeded_hash(data):
    return hmac.new(SEED.encode(), data, hashlib.sha256).hexdigest()


def is_input(path):
    try:
        mode = os.stat(path).st_mode
    except (OSError, ValueError):
        return False
    return stat.S_ISREG(mode) or stat.S_ISFIFO(mode)


def find_inputs(tokens, exclude=(), prefixes=False):
    """Return the input files named by tokens, in order, without duplicates."""
    paths = []
    for token in tokens:
        candidates = [token]
        if token.startswith("-") and "=" in token:
            candidates.append(token.split("=", 1)[1])
        for c in candidates:
            if c in exclude:
                continue
            if is_input(c):
                paths.append(c)
            elif prefixes and c != "":
                import glob
                # Index prefixes like the one used by bwa
                paths.extend(p for p in sorted(glob.glob(glob.escape(c) + ".*")) if is_input(p))
    result = []
    for p in paths:
        if p not in result:
            result.append(p)
    return result


def option(args, name):
    """Return the value of option name, given as "name value" or "name=value"."""
    for i, a in enumerate(args):
        if a == name:
            if i + 1 >= len(args):
                fail(f"option {name} has no value")
            return args[i + 1]
        if a.startswith(name + "="):
            return a[len(name) + 1:]
    fail(f"option {name} is missing")


def record(inputs, stdin=None, extra=()):
    lines = [f"tool\t{TOOL}"]
    lines += [f"arg\t{a}" for a in ARGS]
    lines.append(f"seed\t{SEED}")
    for path in inputs:
        with open(path, "rb") as f:
            lines.append(f"input\t{path}\t{seeded_hash(decompress(f.read()))}")
    if stdin is not None:
        lines.append(f"stdin\t{seeded_hash(stdin)}")
        lines += ["| " + l for l in stdin.decode().splitlines()]
    lines += list(extra)
    return "\n".join(lines) + "\n"


def write(path, text):
    data = text.encode() if isinstance(text, str) else text
    if path.endswith(".gz"):
        import gzip
        data = gzip.compress(data)
    elif path.endswith(".bz2"):
        import bz2
        data = bz2.compress(data)
    with open(path, "wb") as f:
        f.write(data)
    print(f"MOCK {TOOL}: wrote {path}", file=sys.stderr)


def fastq_names(acc):
    layout = os.environ.get("MOCK_SRA_LAYOUT", "paired")
    if layout == "paired":
        return [f"{acc}_1.fastq", f"{acc}_2.fastq"]
    if layout == "single":
        return [f"{acc}.fastq"]
    if layout == "both":
        return [f"{acc}.fastq", f"{acc}_1.fastq", f"{acc}_2.fastq"]
    fail(f"invalid MOCK_SRA_LAYOUT {layout}")


# Tools

def mock_bwa():
    if ARGS[:1] == ["index"]:
        prefix = option(ARGS, "-p")
        write(prefix + ".mockidx", record(find_inputs(ARGS[1:])))
    elif ARGS[:1] == ["mem"]:
        sys.stdout.write(record(find_inputs(ARGS[1:], prefixes=True)))
    else:
        fail(f"unsupported subcommand {ARGS[:1]}")


def mock_samtools():
    if ARGS[:1] == ["sort"]:
        output = option(ARGS, "-o")
        stdin = sys.stdin.buffer.read() if "-" in ARGS[1:] else None
        write(output, record(find_inputs(ARGS[1:], exclude=[output]), stdin=stdin))
    elif ARGS[:1] == ["index"]:
        output = ARGS[-1]
        write(output, record(find_inputs(ARGS[1:-1])))
    else:
        fail(f"unsupported subcommand {ARGS[:1]}")


def mock_fastqc():
    inputs = find_inputs(ARGS)
    if not inputs:
        fail("no input files")
    text = record(inputs)
    for path in inputs:
        base = re.sub(r"(\.gz|\.bz2)$", "", path)
        base = re.sub(r"\.(bam|sam|fastq|fq)$", "", base)
        write(base + "_fastqc.html", text)
        write(base + "_fastqc.zip", text)


def mock_mapping_statistics():
    output = ARGS[-1]
    # The steps read these values from the stats file with grep
    write(output, record(find_inputs(ARGS[:-1]), extra=["insert_size\t250", "read_length\t100"]))


def mock_output_option(name, prefixes=False):
    def mock():
        output = option(ARGS, name)
        write(output, record(find_inputs(ARGS, exclude=[output], prefixes=prefixes)))
    return mock


def mock_gnuplot():
    script = sys.stdin.buffer.read()
    text = script.decode()
    outputs = re.findall(r'set output "([^"]+)"', text)
    if len(outputs) != 1:
        fail(f"expected one 'set output', found {len(outputs)}")
    inputs = find_inputs(re.findall(r'<\s*gzip -dc\s+([^"]+)"', text))
    write(outputs[0], record(inputs, stdin=script))


def mock_rscript():
    # Rscript fit_velocity_oscillation.R <genome> <model> <lambda> <abundance> <output> <cores>
    if len(ARGS) < 6:
        fail("expected at least 6 arguments")
    output = ARGS[5]
    write(output, record(find_inputs(ARGS, exclude=[output])))


def mock_velogp():
    outputs = []
    for i, a in enumerate(ARGS):
        name = a.split("=", 1)[0]
        if name.startswith("--out-") or name.startswith("--plot-"):
            outputs.append(option(ARGS[i:], name))
    if not outputs:
        fail("no outputs")
    text = record(find_inputs(ARGS, exclude=outputs))
    for output in outputs:
        write(output, text)


def mock_featurecounts():
    output = option(ARGS, "-o")
    text = record(find_inputs(ARGS, exclude=[output]))
    write(output, text)
    write(output + ".summary", text)


def mock_prefetch():
    output = option(ARGS, "-o")
    write(output, record(find_inputs(ARGS, exclude=[output])))


def mock_fastq_dump():
    outdir = option(ARGS, "--outdir")
    sra_file = ARGS[-1]
    acc = re.sub(r"\.sra$", "", os.path.basename(sra_file))
    text = record(find_inputs(ARGS))
    for name in fastq_names(acc):
        write(os.path.join(outdir, name), text)


def mock_fasterq_dump():
    acc = ARGS[-1]
    text = record(find_inputs(ARGS))
    for name in fastq_names(acc):
        write(name, text)


def mock_curl():
    output = option(ARGS, "-o")
    urls = [a for a in ARGS if a.startswith("https://")]
    if len(urls) != 1:
        fail(f"expected one URL, found {len(urls)}")
    m = re.search(r"/genome/accession/([^/]+)/download\?", urls[0])
    if m is None:
        fail(f"unexpected URL {urls[0]}")
    acc = m[1]
    types = re.findall(r"include_annotation_type=([A-Z_]+)", urls[0])
    suffixes = {"GENOME_FASTA": f"{acc}_ASM_genomic.fna", "GENOME_GFF": "genomic.gff",
                "GENOME_GBFF": "genomic.gbff"}
    text = record(find_inputs(ARGS, exclude=[output]))
    import io
    import zipfile
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        z.writestr("README.md", "mock NCBI dataset\n")
        for t in types:
            if t not in suffixes:
                fail(f"unsupported annotation type {t}")
            z.writestr(f"ncbi_dataset/data/{acc}/{suffixes[t]}", text)
    write(output, buf.getvalue())


TOOLS = {
    "bwa": mock_bwa,
    "samtools": mock_samtools,
    "fastqc": mock_fastqc,
    "mapping_statistics.py": mock_mapping_statistics,
    "mappability.py": mock_output_option("--output", prefixes=True),
    "abundance.py": mock_output_option("--output"),
    "gnuplot": mock_gnuplot,
    "Rscript": mock_rscript,
    "velogp": mock_velogp,
    "velogp.velomix": mock_velogp,
    "featureCounts": mock_featurecounts,
    "prefetch": mock_prefetch,
    "fastq-dump": mock_fastq_dump,
    "fasterq-dump": mock_fasterq_dump,
    "curl": mock_curl,
}

TOOL = None
ARGS = None
SEED = None


def main(tool=None):
    global TOOL, ARGS, SEED
    TOOL = tool or os.path.basename(sys.argv[0])
    ARGS = sys.argv[1:]
    SEED = os.environ.get("MOCK_SEED")
    if not SEED:
        fail("MOCK_SEED is not set")
    if TOOL not in TOOLS:
        fail("unknown tool")
    TOOLS[TOOL]()


if __name__ == "__main__":
    main()
