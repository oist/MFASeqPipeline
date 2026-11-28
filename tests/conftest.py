"""Helpers for the tests that run pipeline steps with mock tools, see tests/mock/mocktool.py."""
import bz2
import fnmatch
import gzip
import hashlib
import hmac
import os
import random
import subprocess
import sys

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ARUBAITO = os.path.join(REPO, "arubaito.py")
MOCK_BIN = os.path.join(REPO, "tests", "mock", "bin")
MOCK_PIPELINE = os.path.join(REPO, "tests", "mock", "pipeline")

# Options that take no value, per tool. Used to split command lines into options
# and positional arguments; all other options take the argument following them.
FLAGS = {
    "fastqc": {"--noextract"},
    "mappability.py": {"--circular", "--paired"},
    "abundance.py": {"--paired"},
    "velogp": {"--error-model"},
    "featureCounts": {"-p", "--countReadPairs", "-B", "-C", "-P", "-O", "-f"},
    "prefetch": {"--progress"},
    "fastq-dump": {"--split-3", "-v"},
    "fasterq-dump": {"--split-3"},
    "curl": {"--fail"},
}


@pytest.fixture
def seed():
    return str(random.randrange(1, 2**31))


@pytest.fixture
def mock_env(seed):
    env = dict(os.environ)
    # The mock tools come first. The directory of the running Python comes next,
    # so that steps find python3 and python (and the mocks a Python 3).
    env["PATH"] = os.pathsep.join([MOCK_BIN, os.path.dirname(sys.executable), env["PATH"]])
    env["MOCK_SEED"] = seed
    env.pop("MOCK_SRA_LAYOUT", None)
    return env


# Input files

def _write(path, data):
    path = str(path)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if path.endswith(".gz"):
        data = gzip.compress(data)
    elif path.endswith(".bz2"):
        data = bz2.compress(data)
    with open(path, "wb") as f:
        f.write(data)


def write_seed(path, seed):
    """Write an input file that contains just the seed (compressed if the name says so)."""
    _write(path, f"{seed}\n".encode())


def write_seed_dir(path, seed, name="index.seed"):
    """Create an input directory containing one seed file."""
    write_seed(os.path.join(str(path), name), seed)


def write_stats(path, seed, insert_size="250", read_length="100"):
    """Write a stats file. Steps extract insert_size and read_length from it with grep."""
    _write(path, f"{seed}\ninsert_size\t{insert_size}\nread_length\t{read_length}\n".encode())


def write_growthrates(path, seed, sample, rate):
    """Write a growth rates file. Steps extract the rate with awk."""
    _write(path, f"{seed}\n{sample}\t{rate}\n".encode())


# Running arubaito

def run_arubaito(*args, env, cwd):
    """Run arubaito.py with the given arguments, return the CompletedProcess."""
    r = subprocess.run([sys.executable, ARUBAITO] + [str(a) for a in args], env=env, cwd=str(cwd),
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    print(r.stdout)
    return r


# Checking outputs

def read_text(path):
    with open(str(path), "rb") as f:
        data = f.read()
    if data.startswith(b"\x1f\x8b"):
        data = gzip.decompress(data)
    elif data.startswith(b"BZh"):
        data = bz2.decompress(data)
    return data.decode()


def seeded_hash(seed, path):
    with open(str(path), "rb") as f:
        data = f.read()
    if data.startswith(b"\x1f\x8b"):
        data = gzip.decompress(data)
    elif data.startswith(b"BZh"):
        data = bz2.decompress(data)
    return hmac.new(str(seed).encode(), data, hashlib.sha256).hexdigest()


def parse_record(text):
    rec = {"tool": None, "args": [], "seed": None, "inputs": [], "stdin_hash": None,
           "stdin": [], "extra": {}}
    for line in text.splitlines():
        if line.startswith("| "):
            rec["stdin"].append(line[2:])
            continue
        key, _, value = line.partition("\t")
        if key == "tool":
            rec["tool"] = value
        elif key == "arg":
            rec["args"].append(value)
        elif key == "seed":
            rec["seed"] = value
        elif key == "input":
            path, _, h = value.rpartition("\t")
            rec["inputs"].append((path, h))
        elif key == "stdin":
            rec["stdin_hash"] = value
        else:
            rec["extra"][key] = value
    rec["stdin"] = "\n".join(rec["stdin"]) + "\n" if rec["stdin"] else ""
    return rec


def split_args(tool, args):
    """Split a command line into positional arguments and (option, value) pairs."""
    flags = FLAGS.get(tool, set())
    positional, options = [], []
    i = 0
    while i < len(args):
        a = args[i]
        if a.startswith("-") and a not in ("-", "--"):
            if "=" in a and a.startswith("--"):
                options.append(tuple(a.split("=", 1)))
            elif a in flags or i + 1 >= len(args) or (args[i + 1].startswith("-") and args[i + 1] != ""):
                options.append((a,))
            else:
                options.append((a, args[i + 1]))
                i += 1
        else:
            positional.append(a)
        i += 1
    return positional, options


def check_args(tool, actual, expected):
    """Compare command lines. Positional arguments must match in order, options in
    any order. Expected values may contain * wildcards (e.g. for temporary paths)."""
    a_pos, a_opt = split_args(tool, actual)
    e_pos, e_opt = split_args(tool, expected)
    assert len(a_pos) == len(e_pos) and all(fnmatch.fnmatchcase(a, e) for a, e in zip(a_pos, e_pos)), \
        f"positional arguments of {tool} differ:\n  actual:   {a_pos}\n  expected: {e_pos}"
    unmatched = list(a_opt)
    for e in e_opt:
        match = next((a for a in unmatched
                      if len(a) == len(e) and a[0] == e[0] and (len(e) == 1 or fnmatch.fnmatchcase(a[1], e[1]))),
                     None)
        assert match is not None, f"option {e} of {tool} missing, actual options: {a_opt}"
        unmatched.remove(match)
    assert not unmatched, f"unexpected options of {tool}: {unmatched}"


def check_record(text, seed, tool, args, inputs=(), stdin=None, stdin_contains=(), extra=None):
    """Check a record written by a mock tool.

    inputs lists the input files whose hashes the record must contain. An entry is
    either the path of the file, which the tool read under that path or (if inputs
    were copied to the temporary directory) under the same name elsewhere; or a pair
    (recorded path pattern, file) for files the tool read under a different name,
    e.g. through a named pipe. The hash is checked against the content of the file.

    stdin is a dictionary of arguments for check_record() to check the record the
    tool read from its standard input. stdin_contains lists strings that the
    standard input must contain.
    """
    rec = parse_record(text)
    assert rec["tool"] == tool, f"expected a record of {tool}, found:\n{text}"
    assert rec["seed"] == str(seed), f"wrong seed in record of {tool}"
    check_args(tool, rec["args"], args)
    recorded = list(rec["inputs"])
    for entry in inputs:
        if isinstance(entry, tuple):
            pattern, path = str(entry[0]), str(entry[1])
            match = next((r for r in recorded if fnmatch.fnmatchcase(r[0], pattern)), None)
        else:
            path = str(entry)
            match = next((r for r in recorded if r[0] == path), None) or \
                next((r for r in recorded if os.path.basename(r[0]) == os.path.basename(path)), None)
        assert match is not None, f"input {entry} missing in record of {tool}, recorded inputs: {recorded}"
        assert match[1] == seeded_hash(seed, path), f"hash of input {match[0]} of {tool} does not match {path}"
        recorded.remove(match)
    assert not recorded, f"unexpected inputs of {tool}: {recorded}"
    if stdin is not None:
        assert rec["stdin"], f"{tool} did not read its standard input"
        check_record(rec["stdin"], seed, **stdin)
    for s in stdin_contains:
        assert s in rec["stdin"], f"standard input of {tool} does not contain {s!r}:\n{rec['stdin']}"
    for key, value in (extra or {}).items():
        assert rec["extra"].get(key) == value
    return rec


def check_output(path, seed, tool, args, **kwargs):
    """Check that the output file exists and is a record as described for check_record()."""
    assert os.path.isfile(str(path)), f"output {path} missing"
    return check_record(read_text(path), seed, tool, args, **kwargs)
