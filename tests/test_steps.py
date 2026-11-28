"""Tests that run a single pipeline step with arubaito, using the mock tools.

Each test writes the step's input files (containing just the seed), runs the one
step with arubaito directly on a pipeline file, with all settings given on the
command line, and checks that each output was written by the expected tool, with
the expected arguments, from the expected inputs.

Test names must not contain "insert_size" or "read_length": the steps grep for
these in files whose records contain paths below the test's tmp_path.
"""
import os

from conftest import (MOCK_PIPELINE, check_output, check_record, read_text, run_arubaito,
                      write_growthrates, write_seed, write_seed_dir, write_stats)


def mfa_settings(tmp_path, *extra):
    """Settings needed to run any step of mfa_pipeline (and its variants) locally, with mock tools."""
    return [f"mfa_pipeline_dir={MOCK_PIPELINE}", f"base_dir={tmp_path}/data", "project=P", "sample=S",
            f"tmp_base_dir={tmp_path}/tmp", "execute_cmd='{script_file}'", "execute_copy_outputs=False",
            "cores=4", *extra]


def other_settings(tmp_path, *extra):
    """Settings needed to run steps of the pipelines that are used together with mfa_pipeline."""
    return [f"project_dir={tmp_path}/data/P", "sample=S", f"tmp_dir={tmp_path}/tmp/S-", "shell=/bin/bash",
            "execute_cmd='{script_file}'", f"log_file={tmp_path}/data/P/results/S.{{step}}.log",
            "step_cores=4", *extra]


# ---------------------------------------------------------------------------
# genome_index

def test_genome_index(tmp_path, seed, mock_env):
    d = tmp_path / "data"
    write_seed(d / "genome.fasta", seed)

    r = run_arubaito("--steps", "genome_index", "mfa_pipeline", *mfa_settings(tmp_path),
                     env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    check_output(d / "genome.fasta.bwaidx/index.mockidx", seed, "bwa",
                 ["index", "-p", f"{d}/genome.fasta.bwaidx/index", f"{d}/genome.fasta"],
                 inputs=[d / "genome.fasta"])


def test_genome_index_copy_outputs(tmp_path, seed, mock_env):
    d = tmp_path / "data"
    tmp = f"{tmp_path}/tmp/P/S-*"
    write_seed(d / "genome.fasta", seed)

    r = run_arubaito("--steps", "genome_index", "mfa_pipeline", *mfa_settings(tmp_path),
                     "execute_copy_outputs=True", env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    check_output(d / "genome.fasta.bwaidx/index.mockidx", seed, "bwa",
                 ["index", "-p", f"{tmp}/genome_index_dir/genome.fasta.bwaidx/index", f"{d}/genome.fasta"],
                 inputs=[d / "genome.fasta"])


# ---------------------------------------------------------------------------
# map (paired-end)

def test_map_paired_gz(tmp_path, seed, mock_env):
    d, fq, out = tmp_path / "data", tmp_path / "data/P/fastq", tmp_path / "data/P/results"
    tmp = f"{tmp_path}/tmp/P/S-*"
    write_seed(fq / "S_r1.fastq.gz", seed)
    write_seed(fq / "S_r2.fastq.gz", seed)
    write_seed_dir(d / "genome.fasta.bwaidx", seed)

    r = run_arubaito("--steps", "map", "mfa_pipeline", *mfa_settings(tmp_path), env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    check_output(out / "S.bam", seed, "samtools",
                 ["sort", "--threads", "4", "-m", "4096M", "-T", f"{tmp}/S.samtools-sort.", "-O", "BAM",
                  "-l", "7", "-o", f"{out}/S.bam.in-progress", "-"],
                 stdin=dict(tool="bwa",
                            args=["mem", "-t", "4", "-v", "1", f"{d}/genome.fasta.bwaidx/index",
                                  f"{tmp}/S.r1.pipe", f"{tmp}/S.r2.pipe"],
                            inputs=[d / "genome.fasta.bwaidx/index.seed",
                                    (f"{tmp}/S.r1.pipe", fq / "S_r1.fastq.gz"),
                                    (f"{tmp}/S.r2.pipe", fq / "S_r2.fastq.gz")]))
    check_output(out / "S.bai", seed, "samtools",
                 ["index", "-@", "4", f"{out}/S.bam", f"{out}/S.bai"],
                 inputs=[out / "S.bam"])
    assert not os.path.exists(f"{out}/S.bam.in-progress")


def test_map_paired_bz2(tmp_path, seed, mock_env):
    d, fq, out = tmp_path / "data", tmp_path / "data/P/fastq", tmp_path / "data/P/results"
    tmp = f"{tmp_path}/tmp/P/S-*"
    write_seed(fq / "S_r1.fastq.bz2", seed)
    write_seed(fq / "S_r2.fastq.bz2", seed)
    write_seed_dir(d / "genome.fasta.bwaidx", seed)

    r = run_arubaito("--steps", "map", "mfa_pipeline", *mfa_settings(tmp_path),
                     "read1_file='{fastq_dir}/{sample}_r1.fastq.bz2'",
                     "read2_file='{fastq_dir}/{sample}_r2.fastq.bz2'", env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    check_output(out / "S.bam", seed, "samtools",
                 ["sort", "--threads", "4", "-m", "4096M", "-T", f"{tmp}/S.samtools-sort.", "-O", "BAM",
                  "-l", "7", "-o", f"{out}/S.bam.in-progress", "-"],
                 stdin=dict(tool="bwa",
                            args=["mem", "-t", "4", "-v", "1", f"{d}/genome.fasta.bwaidx/index",
                                  f"{tmp}/S.r1.pipe", f"{tmp}/S.r2.pipe"],
                            inputs=[d / "genome.fasta.bwaidx/index.seed",
                                    (f"{tmp}/S.r1.pipe", fq / "S_r1.fastq.bz2"),
                                    (f"{tmp}/S.r2.pipe", fq / "S_r2.fastq.bz2")]))
    check_output(out / "S.bai", seed, "samtools",
                 ["index", "-@", "4", f"{out}/S.bam", f"{out}/S.bai"],
                 inputs=[out / "S.bam"])


def test_map_paired_uncompressed(tmp_path, seed, mock_env):
    d, fq, out = tmp_path / "data", tmp_path / "data/P/fastq", tmp_path / "data/P/results"
    tmp = f"{tmp_path}/tmp/P/S-*"
    write_seed(fq / "S_r1.fastq", seed)
    write_seed(fq / "S_r2.fastq", seed)
    write_seed_dir(d / "genome.fasta.bwaidx", seed)

    r = run_arubaito("--steps", "map", "mfa_pipeline", *mfa_settings(tmp_path),
                     "read1_file='{fastq_dir}/{sample}_r1.fastq'",
                     "read2_file='{fastq_dir}/{sample}_r2.fastq'", env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    check_output(out / "S.bam", seed, "samtools",
                 ["sort", "--threads", "4", "-m", "4096M", "-T", f"{tmp}/S.samtools-sort.", "-O", "BAM",
                  "-l", "7", "-o", f"{out}/S.bam.in-progress", "-"],
                 stdin=dict(tool="bwa",
                            args=["mem", "-t", "4", "-v", "1", f"{d}/genome.fasta.bwaidx/index",
                                  f"{tmp}/S.r1.pipe", f"{tmp}/S.r2.pipe"],
                            inputs=[d / "genome.fasta.bwaidx/index.seed",
                                    (f"{tmp}/S.r1.pipe", fq / "S_r1.fastq"),
                                    (f"{tmp}/S.r2.pipe", fq / "S_r2.fastq")]))
    check_output(out / "S.bai", seed, "samtools",
                 ["index", "-@", "4", f"{out}/S.bam", f"{out}/S.bai"],
                 inputs=[out / "S.bam"])


def test_map_paired_memory_and_cores(tmp_path, seed, mock_env):
    d, fq, out = tmp_path / "data", tmp_path / "data/P/fastq", tmp_path / "data/P/results"
    tmp = f"{tmp_path}/tmp/P/S-*"
    write_seed(fq / "S_r1.fastq.gz", seed)
    write_seed(fq / "S_r2.fastq.gz", seed)
    write_seed_dir(d / "genome.fasta.bwaidx", seed)

    r = run_arubaito("--steps", "map", "mfa_pipeline", *mfa_settings(tmp_path), "cores=8",
                     "memory_sort_gb=2", env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    check_output(out / "S.bam", seed, "samtools",
                 ["sort", "--threads", "8", "-m", "256M", "-T", f"{tmp}/S.samtools-sort.", "-O", "BAM",
                  "-l", "7", "-o", f"{out}/S.bam.in-progress", "-"],
                 stdin=dict(tool="bwa",
                            args=["mem", "-t", "8", "-v", "1", f"{d}/genome.fasta.bwaidx/index",
                                  f"{tmp}/S.r1.pipe", f"{tmp}/S.r2.pipe"],
                            inputs=[d / "genome.fasta.bwaidx/index.seed",
                                    (f"{tmp}/S.r1.pipe", fq / "S_r1.fastq.gz"),
                                    (f"{tmp}/S.r2.pipe", fq / "S_r2.fastq.gz")]))
    check_output(out / "S.bai", seed, "samtools",
                 ["index", "-@", "8", f"{out}/S.bam", f"{out}/S.bai"],
                 inputs=[out / "S.bam"])


def test_map_paired_copy_outputs(tmp_path, seed, mock_env):
    d, fq, out = tmp_path / "data", tmp_path / "data/P/fastq", tmp_path / "data/P/results"
    tmp = f"{tmp_path}/tmp/P/S-*"
    write_seed(fq / "S_r1.fastq.gz", seed)
    write_seed(fq / "S_r2.fastq.gz", seed)
    write_seed_dir(d / "genome.fasta.bwaidx", seed)

    r = run_arubaito("--steps", "map", "mfa_pipeline", *mfa_settings(tmp_path), "execute_copy_outputs=True",
                     env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    check_output(out / "S.bam", seed, "samtools",
                 ["sort", "--threads", "4", "-m", "4096M", "-T", f"{tmp}/S.samtools-sort.", "-O", "BAM",
                  "-l", "7", "-o", f"{tmp}/bam_file/S.bam.in-progress", "-"],
                 stdin=dict(tool="bwa",
                            args=["mem", "-t", "4", "-v", "1", f"{d}/genome.fasta.bwaidx/index",
                                  f"{tmp}/S.r1.pipe", f"{tmp}/S.r2.pipe"],
                            inputs=[d / "genome.fasta.bwaidx/index.seed",
                                    (f"{tmp}/S.r1.pipe", fq / "S_r1.fastq.gz"),
                                    (f"{tmp}/S.r2.pipe", fq / "S_r2.fastq.gz")]))
    check_output(out / "S.bai", seed, "samtools",
                 ["index", "-@", "4", f"{tmp}/bam_file/S.bam", f"{tmp}/bam_file/S.bai"],
                 inputs=[out / "S.bam"])


def test_map_paired_copy_inputs(tmp_path, seed, mock_env):
    d, fq, out = tmp_path / "data", tmp_path / "data/P/fastq", tmp_path / "data/P/results"
    tmp = f"{tmp_path}/tmp/P/S-*"
    write_seed(fq / "S_r1.fastq.gz", seed)
    write_seed(fq / "S_r2.fastq.gz", seed)
    write_seed_dir(d / "genome.fasta.bwaidx", seed)

    r = run_arubaito("--steps", "map", "mfa_pipeline", *mfa_settings(tmp_path), "execute_copy_inputs=True",
                     env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    check_output(out / "S.bam", seed, "samtools",
                 ["sort", "--threads", "4", "-m", "4096M", "-T", f"{tmp}/S.samtools-sort.", "-O", "BAM",
                  "-l", "7", "-o", f"{out}/S.bam.in-progress", "-"],
                 stdin=dict(tool="bwa",
                            args=["mem", "-t", "4", "-v", "1",
                                  f"{tmp}/genome_index_dir/genome.fasta.bwaidx/index",
                                  f"{tmp}/S.r1.pipe", f"{tmp}/S.r2.pipe"],
                            inputs=[d / "genome.fasta.bwaidx/index.seed",
                                    (f"{tmp}/S.r1.pipe", fq / "S_r1.fastq.gz"),
                                    (f"{tmp}/S.r2.pipe", fq / "S_r2.fastq.gz")]))
    check_output(out / "S.bai", seed, "samtools",
                 ["index", "-@", "4", f"{out}/S.bam", f"{out}/S.bai"],
                 inputs=[out / "S.bam"])


# ---------------------------------------------------------------------------
# map (single-end)

def test_map_single_gz(tmp_path, seed, mock_env):
    d, fq, out = tmp_path / "data", tmp_path / "data/P/fastq", tmp_path / "data/P/results"
    tmp = f"{tmp_path}/tmp/P/S-*"
    write_seed(fq / "S.fastq.gz", seed)
    write_seed_dir(d / "genome.fasta.bwaidx", seed)

    r = run_arubaito("--steps", "map", "mfa_se_pipeline", *mfa_settings(tmp_path), env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    check_output(out / "S.bam", seed, "samtools",
                 ["sort", "--threads", "4", "-m", "4096M", "-T", f"{tmp}/S.samtools-sort.", "-O", "BAM",
                  "-l", "7", "-o", f"{out}/S.bam.in-progress", "-"],
                 stdin=dict(tool="bwa",
                            args=["mem", "-t", "4", "-v", "1", f"{d}/genome.fasta.bwaidx/index",
                                  f"{tmp}/S.r.pipe"],
                            inputs=[d / "genome.fasta.bwaidx/index.seed",
                                    (f"{tmp}/S.r.pipe", fq / "S.fastq.gz")]))
    check_output(out / "S.bai", seed, "samtools",
                 ["index", "-@", "4", f"{out}/S.bam", f"{out}/S.bai"],
                 inputs=[out / "S.bam"])


def test_map_single_bz2(tmp_path, seed, mock_env):
    d, fq, out = tmp_path / "data", tmp_path / "data/P/fastq", tmp_path / "data/P/results"
    tmp = f"{tmp_path}/tmp/P/S-*"
    write_seed(fq / "S.fastq.bz2", seed)
    write_seed_dir(d / "genome.fasta.bwaidx", seed)

    r = run_arubaito("--steps", "map", "mfa_se_pipeline", *mfa_settings(tmp_path),
                     "read_file='{fastq_dir}/{sample}.fastq.bz2'", env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    check_output(out / "S.bam", seed, "samtools",
                 ["sort", "--threads", "4", "-m", "4096M", "-T", f"{tmp}/S.samtools-sort.", "-O", "BAM",
                  "-l", "7", "-o", f"{out}/S.bam.in-progress", "-"],
                 stdin=dict(tool="bwa",
                            args=["mem", "-t", "4", "-v", "1", f"{d}/genome.fasta.bwaidx/index",
                                  f"{tmp}/S.r.pipe"],
                            inputs=[d / "genome.fasta.bwaidx/index.seed",
                                    (f"{tmp}/S.r.pipe", fq / "S.fastq.bz2")]))
    check_output(out / "S.bai", seed, "samtools",
                 ["index", "-@", "4", f"{out}/S.bam", f"{out}/S.bai"],
                 inputs=[out / "S.bam"])


def test_map_single_uncompressed(tmp_path, seed, mock_env):
    d, fq, out = tmp_path / "data", tmp_path / "data/P/fastq", tmp_path / "data/P/results"
    tmp = f"{tmp_path}/tmp/P/S-*"
    write_seed(fq / "S.fastq", seed)
    write_seed_dir(d / "genome.fasta.bwaidx", seed)

    r = run_arubaito("--steps", "map", "mfa_se_pipeline", *mfa_settings(tmp_path),
                     "read_file='{fastq_dir}/{sample}.fastq'", env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    check_output(out / "S.bam", seed, "samtools",
                 ["sort", "--threads", "4", "-m", "4096M", "-T", f"{tmp}/S.samtools-sort.", "-O", "BAM",
                  "-l", "7", "-o", f"{out}/S.bam.in-progress", "-"],
                 stdin=dict(tool="bwa",
                            args=["mem", "-t", "4", "-v", "1", f"{d}/genome.fasta.bwaidx/index",
                                  f"{tmp}/S.r.pipe"],
                            inputs=[d / "genome.fasta.bwaidx/index.seed",
                                    (f"{tmp}/S.r.pipe", fq / "S.fastq")]))
    check_output(out / "S.bai", seed, "samtools",
                 ["index", "-@", "4", f"{out}/S.bam", f"{out}/S.bai"],
                 inputs=[out / "S.bam"])


# ---------------------------------------------------------------------------
# fastqc

def test_fastqc(tmp_path, seed, mock_env):
    out = tmp_path / "data/P/results"
    write_seed(out / "S.bam", seed)

    r = run_arubaito("--steps", "fastqc", "mfa_pipeline", *mfa_settings(tmp_path), env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    for path in [out / "S_fastqc.html", out / "S_fastqc.zip"]:
        check_output(path, seed, "fastqc", ["--threads", "4", "--noextract", f"{out}/S.bam"],
                     inputs=[out / "S.bam"])


def test_fastqc_custom_names(tmp_path, seed, mock_env):
    out = tmp_path / "data/P/results"
    write_seed(out / "S.bam", seed)

    r = run_arubaito("--steps", "fastqc", "mfa_pipeline", *mfa_settings(tmp_path),
                     "fastqc_html_file='{out_dir}/qc/{sample}.html'", "fastqc_zip_file='{out_dir}/qc/{sample}.zip'",
                     env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    for path in [out / "qc/S.html", out / "qc/S.zip"]:
        check_output(path, seed, "fastqc", ["--threads", "4", "--noextract", f"{out}/S.bam"],
                     inputs=[out / "S.bam"])
    assert not os.path.exists(out / "S_fastqc.html")
    assert not os.path.exists(out / "S_fastqc.zip")


def test_fastqc_copy_outputs(tmp_path, seed, mock_env):
    out = tmp_path / "data/P/results"
    write_seed(out / "S.bam", seed)

    r = run_arubaito("--steps", "fastqc", "mfa_pipeline", *mfa_settings(tmp_path), "execute_copy_outputs=True",
                     env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    for path in [out / "S_fastqc.html", out / "S_fastqc.zip"]:
        check_output(path, seed, "fastqc", ["--threads", "4", "--noextract", f"{out}/S.bam"],
                     inputs=[out / "S.bam"])


# ---------------------------------------------------------------------------
# mapstats

def test_mapstats(tmp_path, seed, mock_env):
    out = tmp_path / "data/P/results"
    write_seed(out / "S.bam", seed)
    write_seed(out / "S.bai", seed)

    r = run_arubaito("--steps", "mapstats", "mfa_pipeline", *mfa_settings(tmp_path), env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    check_output(out / "S.stats", seed, "mapping_statistics.py", [f"{out}/S.bam", f"{out}/S.stats"],
                 inputs=[out / "S.bam"])


def test_mapstats_single(tmp_path, seed, mock_env):
    out = tmp_path / "data/P/results"
    write_seed(out / "S.bam", seed)
    write_seed(out / "S.bai", seed)

    r = run_arubaito("--steps", "mapstats", "mfa_se_pipeline", *mfa_settings(tmp_path), env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    check_output(out / "S.stats", seed, "mapping_statistics.py", [f"{out}/S.bam", f"{out}/S.stats"],
                 inputs=[out / "S.bam"])


# ---------------------------------------------------------------------------
# mappability (paired-end)

def test_mappability_circular_paired(tmp_path, seed, mock_env):
    d, out = tmp_path / "data", tmp_path / "data/P/results"
    write_seed(d / "genome.fasta", seed)
    write_seed_dir(d / "genome.fasta.bwaidx", seed)
    write_stats(out / "S.stats", seed, insert_size="300", read_length="150")

    r = run_arubaito("--steps", "mappability", "mfa_pipeline", *mfa_settings(tmp_path), "genome_circular=True",
                     env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    check_output(out / "S.mappability.bed", seed, "mappability.py",
                 ["--genome", f"{d}/genome.fasta", "--genome-index", f"{d}/genome.fasta.bwaidx/index",
                  "--circular", "--paired", "--fragment-length", "300", "--read-length", "150",
                  "--mapping-quality", "20", "--regions", "", "--exclude-regions", "",
                  "--output", f"{out}/S.mappability.bed", "--", "-t", "4", "-v", "1"],
                 inputs=[d / "genome.fasta", d / "genome.fasta.bwaidx/index.seed"])


def test_mappability_linear_paired(tmp_path, seed, mock_env):
    d, out = tmp_path / "data", tmp_path / "data/P/results"
    write_seed(d / "genome.fasta", seed)
    write_seed_dir(d / "genome.fasta.bwaidx", seed)
    write_stats(out / "S.stats", seed, insert_size="300", read_length="150")

    r = run_arubaito("--steps", "mappability", "mfa_pipeline", *mfa_settings(tmp_path), "genome_circular=False",
                     env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    check_output(out / "S.mappability.bed", seed, "mappability.py",
                 ["--genome", f"{d}/genome.fasta", "--genome-index", f"{d}/genome.fasta.bwaidx/index",
                  "--paired", "--fragment-length", "300", "--read-length", "150",
                  "--mapping-quality", "20", "--regions", "", "--exclude-regions", "",
                  "--output", f"{out}/S.mappability.bed", "--", "-t", "4", "-v", "1"],
                 inputs=[d / "genome.fasta", d / "genome.fasta.bwaidx/index.seed"])


def test_mappability_regions_and_quality(tmp_path, seed, mock_env):
    d, out = tmp_path / "data", tmp_path / "data/P/results"
    write_seed(d / "genome.fasta", seed)
    write_seed_dir(d / "genome.fasta.bwaidx", seed)
    write_seed(d / "regions.bed.gz", seed)
    write_stats(out / "S.stats", seed, insert_size="300", read_length="150")

    r = run_arubaito("--steps", "mappability", "mfa_pipeline", *mfa_settings(tmp_path),
                     "regions_file='{base_dir}/regions.bed.gz'", "exclude_regions=rRNA,phage",
                     "mapping_quality=30", env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    check_output(out / "S.mappability.bed", seed, "mappability.py",
                 ["--genome", f"{d}/genome.fasta", "--genome-index", f"{d}/genome.fasta.bwaidx/index",
                  "--circular", "--paired", "--fragment-length", "300", "--read-length", "150",
                  "--mapping-quality", "30", "--regions", f"{d}/regions.bed.gz", "--exclude-regions", "rRNA,phage",
                  "--output", f"{out}/S.mappability.bed", "--", "-t", "4", "-v", "1"],
                 inputs=[d / "genome.fasta", d / "genome.fasta.bwaidx/index.seed", d / "regions.bed.gz"])


def test_mappability_invalid_circular(tmp_path, seed, mock_env):
    d, out = tmp_path / "data", tmp_path / "data/P/results"
    write_seed(d / "genome.fasta", seed)
    write_seed_dir(d / "genome.fasta.bwaidx", seed)
    write_stats(out / "S.stats", seed)

    r = run_arubaito("--steps", "mappability", "mfa_pipeline", *mfa_settings(tmp_path), "genome_circular=yes please",
                     env=mock_env, cwd=tmp_path)
    assert r.returncode != 0

    assert "genome_circular must be True or False, not 'yes please'" in read_text(out / "S.mappability.log")
    assert not os.path.exists(out / "S.mappability.bed")


def test_mappability_invalid_paired(tmp_path, seed, mock_env):
    d, out = tmp_path / "data", tmp_path / "data/P/results"
    write_seed(d / "genome.fasta", seed)
    write_seed_dir(d / "genome.fasta.bwaidx", seed)
    write_stats(out / "S.stats", seed)

    r = run_arubaito("--steps", "mappability", "mfa_pipeline", *mfa_settings(tmp_path), "paired_end=maybe",
                     env=mock_env, cwd=tmp_path)
    assert r.returncode != 0

    assert "paired_end must be True or False, not 'maybe'" in read_text(out / "S.mappability.log")
    assert not os.path.exists(out / "S.mappability.bed")


# ---------------------------------------------------------------------------
# mappability (single-end)

def test_mappability_circular_single(tmp_path, seed, mock_env):
    d, out = tmp_path / "data", tmp_path / "data/P/results"
    write_seed(d / "genome.fasta", seed)
    write_seed_dir(d / "genome.fasta.bwaidx", seed)
    write_stats(out / "S.stats", seed, insert_size="nan", read_length="150")

    r = run_arubaito("--steps", "mappability", "mfa_se_pipeline", *mfa_settings(tmp_path), "genome_circular=True",
                     env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    check_output(out / "S.mappability.bed", seed, "mappability.py",
                 ["--genome", f"{d}/genome.fasta", "--genome-index", f"{d}/genome.fasta.bwaidx/index",
                  "--circular", "--read-length", "150",
                  "--mapping-quality", "20", "--regions", "", "--exclude-regions", "",
                  "--output", f"{out}/S.mappability.bed", "--", "-t", "4", "-v", "1"],
                 inputs=[d / "genome.fasta", d / "genome.fasta.bwaidx/index.seed"])


def test_mappability_linear_single(tmp_path, seed, mock_env):
    d, out = tmp_path / "data", tmp_path / "data/P/results"
    write_seed(d / "genome.fasta", seed)
    write_seed_dir(d / "genome.fasta.bwaidx", seed)
    write_stats(out / "S.stats", seed, insert_size="nan", read_length="150")

    r = run_arubaito("--steps", "mappability", "mfa_se_pipeline", *mfa_settings(tmp_path), "genome_circular=False",
                     env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    check_output(out / "S.mappability.bed", seed, "mappability.py",
                 ["--genome", f"{d}/genome.fasta", "--genome-index", f"{d}/genome.fasta.bwaidx/index",
                  "--read-length", "150",
                  "--mapping-quality", "20", "--regions", "", "--exclude-regions", "",
                  "--output", f"{out}/S.mappability.bed", "--", "-t", "4", "-v", "1"],
                 inputs=[d / "genome.fasta", d / "genome.fasta.bwaidx/index.seed"])


def test_mappability_single_regions(tmp_path, seed, mock_env):
    d, out = tmp_path / "data", tmp_path / "data/P/results"
    write_seed(d / "genome.fasta", seed)
    write_seed_dir(d / "genome.fasta.bwaidx", seed)
    write_seed(d / "regions.bed", seed)
    write_stats(out / "S.stats", seed, insert_size="nan", read_length="150")

    r = run_arubaito("--steps", "mappability", "mfa_se_pipeline", *mfa_settings(tmp_path),
                     "regions_file='{base_dir}/regions.bed'", "exclude_regions=rRNA", env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    check_output(out / "S.mappability.bed", seed, "mappability.py",
                 ["--genome", f"{d}/genome.fasta", "--genome-index", f"{d}/genome.fasta.bwaidx/index",
                  "--circular", "--read-length", "150",
                  "--mapping-quality", "20", "--regions", f"{d}/regions.bed", "--exclude-regions", "rRNA",
                  "--output", f"{out}/S.mappability.bed", "--", "-t", "4", "-v", "1"],
                 inputs=[d / "genome.fasta", d / "genome.fasta.bwaidx/index.seed", d / "regions.bed"])


def test_mappability_single_invalid_circular(tmp_path, seed, mock_env):
    d, out = tmp_path / "data", tmp_path / "data/P/results"
    write_seed(d / "genome.fasta", seed)
    write_seed_dir(d / "genome.fasta.bwaidx", seed)
    write_stats(out / "S.stats", seed)

    r = run_arubaito("--steps", "mappability", "mfa_se_pipeline", *mfa_settings(tmp_path), "genome_circular=1",
                     env=mock_env, cwd=tmp_path)
    assert r.returncode != 0

    assert "genome_circular must be True or False, not '1'" in read_text(out / "S.mappability.log")
    assert not os.path.exists(out / "S.mappability.bed")


# ---------------------------------------------------------------------------
# abundance

def test_abundance_paired(tmp_path, seed, mock_env):
    out = tmp_path / "data/P/results"
    for f in ["S.bam", "S.bai", "S.mappability.bed", "R.bam", "R.bai", "R.mappability.bed"]:
        write_seed(out / f, seed)

    r = run_arubaito("--steps", "abundance", "mfa_pipeline", *mfa_settings(tmp_path), "reference=R",
                     env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    check_output(out / "S.abd.w10kbp.gz", seed, "abundance.py",
                 ["--sample-bam", f"{out}/S.bam", "--sample-mappability", f"{out}/S.mappability.bed",
                  "--reference-bam", f"{out}/R.bam", "--reference-mappability", f"{out}/R.mappability.bed",
                  "--paired", "--mapping-quality", "20", "--windowsize", "10k", "--output", f"{out}/S.abd.w10kbp.gz"],
                 inputs=[out / "S.bam", out / "S.mappability.bed", out / "R.bam", out / "R.mappability.bed"])
    assert os.path.exists(out / "S.abundance.w10kbp.log")


def test_abundance_paired_windowsize_and_quality(tmp_path, seed, mock_env):
    out = tmp_path / "data/P/results"
    for f in ["S.bam", "S.bai", "S.mappability.bed", "R.bam", "R.bai", "R.mappability.bed"]:
        write_seed(out / f, seed)

    r = run_arubaito("--steps", "abundance", "mfa_pipeline", *mfa_settings(tmp_path), "reference=R",
                     "windowsize=500", "mapping_quality=0", env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    check_output(out / "S.abd.w500bp.gz", seed, "abundance.py",
                 ["--sample-bam", f"{out}/S.bam", "--sample-mappability", f"{out}/S.mappability.bed",
                  "--reference-bam", f"{out}/R.bam", "--reference-mappability", f"{out}/R.mappability.bed",
                  "--paired", "--mapping-quality", "0", "--windowsize", "500", "--output", f"{out}/S.abd.w500bp.gz"],
                 inputs=[out / "S.bam", out / "S.mappability.bed", out / "R.bam", out / "R.mappability.bed"])
    assert os.path.exists(out / "S.abundance.w500bp.log")


def test_abundance_paired_copy_inputs_and_outputs(tmp_path, seed, mock_env):
    out = tmp_path / "data/P/results"
    tmp = f"{tmp_path}/tmp/P/S-*"
    for f in ["S.bam", "S.bai", "S.mappability.bed", "R.bam", "R.bai", "R.mappability.bed"]:
        write_seed(out / f, seed)

    r = run_arubaito("--steps", "abundance", "mfa_pipeline", *mfa_settings(tmp_path), "reference=R",
                     "execute_copy_inputs=True", "execute_copy_outputs=True", env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    # All inputs are in the same directory, so arubaito copies them all to the directory of the first
    check_output(out / "S.abd.w10kbp.gz", seed, "abundance.py",
                 ["--sample-bam", f"{tmp}/bam_file/S.bam", "--sample-mappability", f"{tmp}/bam_file/S.mappability.bed",
                  "--reference-bam", f"{tmp}/bam_file/R.bam",
                  "--reference-mappability", f"{tmp}/bam_file/R.mappability.bed",
                  "--paired", "--mapping-quality", "20", "--windowsize", "10k",
                  "--output", f"{tmp}/bam_file/S.abd.w10kbp.gz"],
                 inputs=[out / "S.bam", out / "S.mappability.bed", out / "R.bam", out / "R.mappability.bed"])


def test_abundance_single(tmp_path, seed, mock_env):
    out = tmp_path / "data/P/results"
    for f in ["S.bam", "S.bai", "S.mappability.bed", "R.bam", "R.bai", "R.mappability.bed"]:
        write_seed(out / f, seed)

    r = run_arubaito("--steps", "abundance", "mfa_se_pipeline", *mfa_settings(tmp_path), "reference=R",
                     env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    check_output(out / "S.abd.w10kbp.gz", seed, "abundance.py",
                 ["--sample-bam", f"{out}/S.bam", "--sample-mappability", f"{out}/S.mappability.bed",
                  "--reference-bam", f"{out}/R.bam", "--reference-mappability", f"{out}/R.mappability.bed",
                  "--mapping-quality", "20", "--windowsize", "10k", "--output", f"{out}/S.abd.w10kbp.gz"],
                 inputs=[out / "S.bam", out / "S.mappability.bed", out / "R.bam", out / "R.mappability.bed"])


def test_abundance_no_reference(tmp_path, seed, mock_env):
    out = tmp_path / "data/P/results"
    for f in ["S.bam", "S.bai", "S.mappability.bed"]:
        write_seed(out / f, seed)

    # mfa_no_reference is an overlay that must be used together with mfa_pipeline
    (tmp_path / "noref.yaml").write_text("use: [mfa_no_reference, mfa_pipeline]\n")

    r = run_arubaito("--steps", "abundance", "noref", *mfa_settings(tmp_path),
                     env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    check_output(out / "S.abd.w10kbp.gz", seed, "abundance.py",
                 ["--sample-bam", f"{out}/S.bam", "--sample-mappability", f"{out}/S.mappability.bed",
                  "--paired", "--mapping-quality", "20", "--windowsize", "10k", "--output", f"{out}/S.abd.w10kbp.gz"],
                 inputs=[out / "S.bam", out / "S.mappability.bed"])


def test_abundance_no_reference_single(tmp_path, seed, mock_env):
    out = tmp_path / "data/P/results"
    for f in ["S.bam", "S.bai", "S.mappability.bed"]:
        write_seed(out / f, seed)

    # mfa_no_reference is an overlay that must be used together with mfa_pipeline
    (tmp_path / "noref.yaml").write_text("use: [mfa_no_reference, mfa_pipeline]\n")

    r = run_arubaito("--steps", "abundance", "noref", *mfa_settings(tmp_path), "paired_end=False",
                     env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    check_output(out / "S.abd.w10kbp.gz", seed, "abundance.py",
                 ["--sample-bam", f"{out}/S.bam", "--sample-mappability", f"{out}/S.mappability.bed",
                  "--mapping-quality", "20", "--windowsize", "10k", "--output", f"{out}/S.abd.w10kbp.gz"],
                 inputs=[out / "S.bam", out / "S.mappability.bed"])


def test_abundance_invalid_paired(tmp_path, seed, mock_env):
    out = tmp_path / "data/P/results"
    for f in ["S.bam", "S.bai", "S.mappability.bed", "R.bam", "R.bai", "R.mappability.bed"]:
        write_seed(out / f, seed)

    r = run_arubaito("--steps", "abundance", "mfa_pipeline", *mfa_settings(tmp_path), "reference=R",
                     "paired_end=paired", env=mock_env, cwd=tmp_path)
    assert r.returncode != 0

    assert "paired_end must be True or False, not 'paired'" in read_text(out / "S.abundance.w10kbp.log")
    assert not os.path.exists(out / "S.abd.w10kbp.gz")


# ---------------------------------------------------------------------------
# abundance_plot

def test_abundance_plot(tmp_path, seed, mock_env):
    out = tmp_path / "data/P/results"
    write_seed(out / "S.abd.w10kbp.gz", seed)

    r = run_arubaito("--steps", "abundance_plot", "mfa_pipeline", *mfa_settings(tmp_path),
                     env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    check_output(out / "S.abd.w10kbp.pdf", seed, "gnuplot", [],
                 inputs=[out / "S.abd.w10kbp.gz"],
                 stdin_contains=["set terminal pdf", f'set output "{out}/S.abd.w10kbp.pdf"',
                                 'column("coverage_ratio")', 'column("ratio_mappable")', 'title "S (well mappable'])


def test_abundance_plot_no_reference(tmp_path, seed, mock_env):
    out = tmp_path / "data/P/results"
    write_seed(out / "S.abd.w5kbp.gz", seed)

    # mfa_no_reference is an overlay that must be used together with mfa_pipeline
    (tmp_path / "noref.yaml").write_text("use: [mfa_no_reference, mfa_pipeline]\n")

    r = run_arubaito("--steps", "abundance_plot", "noref", *mfa_settings(tmp_path), "windowsize=5k",
                     env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    rec = check_output(out / "S.abd.w5kbp.pdf", seed, "gnuplot", [],
                       inputs=[out / "S.abd.w5kbp.gz"],
                       stdin_contains=["set terminal pdf", f'set output "{out}/S.abd.w5kbp.pdf"',
                                       'column("sample_coverage_avg")', 'column("sample_mappable")'])
    assert "coverage_ratio" not in rec["stdin"]


# ---------------------------------------------------------------------------
# velocity_oscillation

def test_velocity_oscillation_vt(tmp_path, seed, mock_env):
    out, meta = tmp_path / "data/P/results", tmp_path / "data/P/meta"
    write_seed(out / "S.abd.w10kbp.gz", seed)
    write_growthrates(meta / "growthrates.tab", seed, "S", "0.693")

    r = run_arubaito("--steps", "velocity_oscillation", "mfa_pipeline", *mfa_settings(tmp_path),
                     "genome_name=EColiMG1655", env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    check_output(out / "S.osc.vt.tab", seed, "Rscript",
                 [f"{MOCK_PIPELINE}/steps/fit_velocity_oscillation.R", "EColiMG1655", "vt", "0.693",
                  f"{out}/S.abd.w10kbp.gz", f"{out}/S.osc.vt.tab", "4"],
                 inputs=[out / "S.abd.w10kbp.gz"])


def test_velocity_oscillation_vx(tmp_path, seed, mock_env):
    out, meta = tmp_path / "data/P/results", tmp_path / "data/P/meta"
    write_seed(out / "S.abd.w10kbp.gz", seed)
    write_growthrates(meta / "growthrates.tab", seed, "S", "1.386")

    r = run_arubaito("--steps", "velocity_oscillation", "mfa_pipeline", *mfa_settings(tmp_path),
                     "genome_name=EColiMG1655", "velocity_oscillation_model=vx", env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    check_output(out / "S.osc.vx.tab", seed, "Rscript",
                 [f"{MOCK_PIPELINE}/steps/fit_velocity_oscillation.R", "EColiMG1655", "vx", "1.386",
                  f"{out}/S.abd.w10kbp.gz", f"{out}/S.osc.vx.tab", "4"],
                 inputs=[out / "S.abd.w10kbp.gz"])


def test_velocity_oscillation_vd(tmp_path, seed, mock_env):
    out, meta = tmp_path / "data/P/results", tmp_path / "data/P/meta"
    write_seed(out / "S.abd.w10kbp.gz", seed)
    write_growthrates(meta / "growthrates.tab", seed, "S", "0.5")

    r = run_arubaito("--steps", "velocity_oscillation", "mfa_pipeline", *mfa_settings(tmp_path),
                     "genome_name=EColiMG1655", "velocity_oscillation_model=vd", env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    check_output(out / "S.osc.vd.tab", seed, "Rscript",
                 [f"{MOCK_PIPELINE}/steps/fit_velocity_oscillation.R", "EColiMG1655", "vd", "0.5",
                  f"{out}/S.abd.w10kbp.gz", f"{out}/S.osc.vd.tab", "4"],
                 inputs=[out / "S.abd.w10kbp.gz"])


# ---------------------------------------------------------------------------
# velogp

def test_velogp(tmp_path, seed, mock_env):
    d, out, meta = tmp_path / "data", tmp_path / "data/P/results", tmp_path / "data/P/meta"
    write_seed(out / "S.abd.w10kbp.gz", seed)
    write_seed(d / "genome.fasta", seed)
    write_seed(d / "genome.gff", seed)
    write_growthrates(meta / "growthrates.tab", seed, "S", "0.693")

    r = run_arubaito("--steps", "velogp", "mfa_pipeline", *mfa_settings(tmp_path), env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    for name in ["abd.gz", "abd.pdf", "vt.gz", "vt.pdf", "vx.gz", "pkl"]:
        check_output(out / f"S.velogp.{name}", seed, "velogp",
                     ["--cores", "4", "--genome", f"{d}/genome.fasta", "--gff", f"{d}/genome.gff",
                      "--columns", "seq,midpoint,coverage_ratio,coverage_ratio_stderr,ratio_mappable",
                      "--vx-feature-types", "", "--vx-feature-names", "", "--vx-delayscale", "0.01",
                      "--vt-covariance", "m1", "--vt-timescale", "0.2", "--vt-magnitude", "0.25",
                      "--out-abundance", f"{out}/S.velogp.abd.gz", "--out-trajectories", f"{out}/S.velogp.vt.gz",
                      "--out-vx-feature-delays", f"{out}/S.velogp.vx.gz", "--out-objects", f"{out}/S.velogp.pkl",
                      "--plot-vt", f"{out}/S.velogp.vt.pdf", "--plot-abundance", f"{out}/S.velogp.abd.pdf",
                      "--origins", "origin_of_replication", "--growthrate", "0.693", "--seed", "42",
                      "--perturbation-rejects", "5", "--perturbation-reduction", "0.8408964",
                      "--perturbation-min", "0.01", f"{out}/S.abd.w10kbp.gz"],
                     inputs=[out / "S.abd.w10kbp.gz", d / "genome.fasta", d / "genome.gff"])
    assert os.path.exists(out / "S.velogp.log")


def test_velogp_all_settings(tmp_path, seed, mock_env):
    d, out, meta = tmp_path / "data", tmp_path / "data/P/results", tmp_path / "data/P/meta"
    write_seed(out / "S.abd.w5kbp.gz", seed)
    write_seed(d / "genome.fasta", seed)
    write_seed(d / "genome.gff", seed)
    write_seed(d / "origins.bed", seed)
    write_growthrates(meta / "growthrates.tab", seed, "S", "1.2")

    r = run_arubaito("--steps", "velogp", "mfa_pipeline", *mfa_settings(tmp_path),
                     "windowsize=5k", "origin_file='{base_dir}/origins.bed'", "origins=oriC,chr2:1000",
                     "velogp_vx_types=tRNA,rRNA", "velogp_vx_names=terA", "velogp_vx_delayscale=0.02",
                     "velogp_vt_cov=se", "velogp_vt_timescale=0.3", "velogp_vt_magnitude=0.5",
                     "velogp_suffix=.test", "velogp_seed=7", "velogp_perturbation_rejects=3",
                     "velogp_perturbation_reduction=0.5", "velogp_perturbation_min=0.001",
                     "velogp_errormodel=True", "velogp_extra_args=--vt-tmax 5 --verbose",
                     env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    for name in ["abd.gz", "abd.pdf", "vt.gz", "vt.pdf", "vx.gz", "pkl"]:
        check_output(out / f"S.velogp.test.{name}", seed, "velogp",
                     ["--cores", "4", "--genome", f"{d}/genome.fasta", "--gff", f"{d}/genome.gff",
                      "--columns", "seq,midpoint,coverage_ratio,coverage_ratio_stderr,ratio_mappable",
                      "--vx-feature-types", "tRNA,rRNA", "--vx-feature-names", "terA", "--vx-delayscale", "0.02",
                      "--vt-covariance", "se", "--vt-timescale", "0.3", "--vt-magnitude", "0.5",
                      "--out-abundance", f"{out}/S.velogp.test.abd.gz",
                      "--out-trajectories", f"{out}/S.velogp.test.vt.gz",
                      "--out-vx-feature-delays", f"{out}/S.velogp.test.vx.gz",
                      "--out-objects", f"{out}/S.velogp.test.pkl",
                      "--plot-vt", f"{out}/S.velogp.test.vt.pdf", "--plot-abundance", f"{out}/S.velogp.test.abd.pdf",
                      "--origins", "oriC,chr2:1000", "--growthrate", "1.2", "--seed", "7",
                      "--perturbation-rejects", "3", "--perturbation-reduction", "0.5",
                      "--perturbation-min", "0.001", "--bed", f"{d}/origins.bed", "--error-model",
                      "--vt-tmax", "5", "--verbose", f"{out}/S.abd.w5kbp.gz"],
                     inputs=[out / "S.abd.w5kbp.gz", d / "genome.fasta", d / "genome.gff", d / "origins.bed"])
    assert os.path.exists(out / "S.velogp.test.log")


def test_velogp_no_reference(tmp_path, seed, mock_env):
    d, out, meta = tmp_path / "data", tmp_path / "data/P/results", tmp_path / "data/P/meta"
    write_seed(out / "S.abd.w10kbp.gz", seed)
    write_seed(d / "genome.fasta", seed)
    write_seed(d / "genome.gff", seed)
    write_growthrates(meta / "growthrates.tab", seed, "S", "0.693")

    # mfa_no_reference is an overlay that must be used together with mfa_pipeline
    (tmp_path / "noref.yaml").write_text("use: [mfa_no_reference, mfa_pipeline]\n")

    r = run_arubaito("--steps", "velogp", "noref", *mfa_settings(tmp_path), env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    for name in ["abd.gz", "abd.pdf", "vt.gz", "vt.pdf", "vx.gz", "pkl"]:
        check_output(out / f"S.velogp.{name}", seed, "velogp",
                     ["--cores", "4", "--genome", f"{d}/genome.fasta", "--gff", f"{d}/genome.gff",
                      "--columns", "seq,midpoint,sample_coverage_avg,sample_coverage_avg_stderr,sample_mappable",
                      "--vx-feature-types", "", "--vx-feature-names", "", "--vx-delayscale", "0.01",
                      "--vt-covariance", "m1", "--vt-timescale", "0.2", "--vt-magnitude", "0.25",
                      "--out-abundance", f"{out}/S.velogp.abd.gz", "--out-trajectories", f"{out}/S.velogp.vt.gz",
                      "--out-vx-feature-delays", f"{out}/S.velogp.vx.gz", "--out-objects", f"{out}/S.velogp.pkl",
                      "--plot-vt", f"{out}/S.velogp.vt.pdf", "--plot-abundance", f"{out}/S.velogp.abd.pdf",
                      "--origins", "origin_of_replication", "--growthrate", "0.693", "--seed", "42",
                      "--perturbation-rejects", "5", "--perturbation-reduction", "0.8408964",
                      "--perturbation-min", "0.01", f"{out}/S.abd.w10kbp.gz"],
                     inputs=[out / "S.abd.w10kbp.gz", d / "genome.fasta", d / "genome.gff"])


def test_velogp_invalid_errormodel(tmp_path, seed, mock_env):
    d, out, meta = tmp_path / "data", tmp_path / "data/P/results", tmp_path / "data/P/meta"
    write_seed(out / "S.abd.w10kbp.gz", seed)
    write_seed(d / "genome.fasta", seed)
    write_seed(d / "genome.gff", seed)
    write_growthrates(meta / "growthrates.tab", seed, "S", "0.693")

    r = run_arubaito("--steps", "velogp", "mfa_pipeline", *mfa_settings(tmp_path), "velogp_errormodel=enabled",
                     env=mock_env, cwd=tmp_path)
    assert r.returncode != 0

    assert "velogp_errormodel must be True or False, not 'enabled'" in read_text(out / "S.velogp.log")
    assert not os.path.exists(out / "S.velogp.pkl")


def test_velogp_copy_outputs(tmp_path, seed, mock_env):
    d, out, meta = tmp_path / "data", tmp_path / "data/P/results", tmp_path / "data/P/meta"
    tmp = f"{tmp_path}/tmp/P/S-*"
    write_seed(out / "S.abd.w10kbp.gz", seed)
    write_seed(d / "genome.fasta", seed)
    write_seed(d / "genome.gff", seed)
    write_growthrates(meta / "growthrates.tab", seed, "S", "0.693")

    r = run_arubaito("--steps", "velogp", "mfa_pipeline", *mfa_settings(tmp_path), "execute_copy_outputs=True",
                     env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    o = f"{tmp}/velogp_abundance_file"
    for name in ["abd.gz", "abd.pdf", "vt.gz", "vt.pdf", "vx.gz", "pkl"]:
        check_output(out / f"S.velogp.{name}", seed, "velogp",
                     ["--cores", "4", "--genome", f"{d}/genome.fasta", "--gff", f"{d}/genome.gff",
                      "--columns", "seq,midpoint,coverage_ratio,coverage_ratio_stderr,ratio_mappable",
                      "--vx-feature-types", "", "--vx-feature-names", "", "--vx-delayscale", "0.01",
                      "--vt-covariance", "m1", "--vt-timescale", "0.2", "--vt-magnitude", "0.25",
                      "--out-abundance", f"{o}/S.velogp.abd.gz", "--out-trajectories", f"{o}/S.velogp.vt.gz",
                      "--out-vx-feature-delays", f"{o}/S.velogp.vx.gz", "--out-objects", f"{o}/S.velogp.pkl",
                      "--plot-vt", f"{o}/S.velogp.vt.pdf", "--plot-abundance", f"{o}/S.velogp.abd.pdf",
                      "--origins", "origin_of_replication", "--growthrate", "0.693", "--seed", "42",
                      "--perturbation-rejects", "5", "--perturbation-reduction", "0.8408964",
                      "--perturbation-min", "0.01", f"{out}/S.abd.w10kbp.gz"],
                     inputs=[out / "S.abd.w10kbp.gz", d / "genome.fasta", d / "genome.gff"])


# ---------------------------------------------------------------------------
# velomix

def test_velomix(tmp_path, seed, mock_env):
    out = tmp_path / "data/P/results"
    write_seed(out / "S.velogp.pkl", seed)

    r = run_arubaito("--steps", "velomix", "mfa_pipeline", *mfa_settings(tmp_path), env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    for name in ["abd.gz", "abd.pdf", "p.gz", "vt.gz", "vt.pdf", "vx.gz", "pkl"]:
        check_output(out / f"S.velomix.{name}", seed, "velogp.velomix",
                     ["--stall-sites", "", "--cores", "4",
                      "--out-probabilities", f"{out}/S.velomix.p.gz", "--out-abundance", f"{out}/S.velomix.abd.gz",
                      "--out-trajectories", f"{out}/S.velomix.vt.gz",
                      "--out-vx-feature-delays", f"{out}/S.velomix.vx.gz", "--out-objects", f"{out}/S.velomix.pkl",
                      "--plot-vt", f"{out}/S.velomix.vt.pdf", "--plot-abundance", f"{out}/S.velomix.abd.pdf",
                      "--seed", "42", "--perturbation-rejects", "5", "--perturbation-reduction", "0.8408964",
                      "--perturbation-min", "0.01", f"{out}/S.velogp.pkl"],
                     inputs=[out / "S.velogp.pkl"])
    assert os.path.exists(out / "S.velomix.log")


def test_velomix_all_settings(tmp_path, seed, mock_env):
    out = tmp_path / "data/P/results"
    write_seed(out / "S.velogp.test.pkl", seed)

    r = run_arubaito("--steps", "velomix", "mfa_pipeline", *mfa_settings(tmp_path),
                     "velogp_suffix=.test", "velomix_suffix=.ter", "velomix_stall_sites=terA,terB",
                     "velomix_seed=3", "velomix_perturbation_rejects=2", "velomix_perturbation_reduction=0.7",
                     "velomix_perturbation_min=0.1", "velomix_extra_args=--verbose", env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    for name in ["abd.gz", "abd.pdf", "p.gz", "vt.gz", "vt.pdf", "vx.gz", "pkl"]:
        check_output(out / f"S.velomix.ter.{name}", seed, "velogp.velomix",
                     ["--stall-sites", "terA,terB", "--cores", "4",
                      "--out-probabilities", f"{out}/S.velomix.ter.p.gz",
                      "--out-abundance", f"{out}/S.velomix.ter.abd.gz",
                      "--out-trajectories", f"{out}/S.velomix.ter.vt.gz",
                      "--out-vx-feature-delays", f"{out}/S.velomix.ter.vx.gz",
                      "--out-objects", f"{out}/S.velomix.ter.pkl",
                      "--plot-vt", f"{out}/S.velomix.ter.vt.pdf", "--plot-abundance", f"{out}/S.velomix.ter.abd.pdf",
                      "--seed", "3", "--perturbation-rejects", "2", "--perturbation-reduction", "0.7",
                      "--perturbation-min", "0.1", "--verbose", f"{out}/S.velogp.test.pkl"],
                     inputs=[out / "S.velogp.test.pkl"])
    assert os.path.exists(out / "S.velomix.ter.log")


# ---------------------------------------------------------------------------
# genecounts (rnaseq_pipeline)

def test_genecounts(tmp_path, seed, mock_env):
    d, out = tmp_path / "data", tmp_path / "data/P/results"
    write_seed(out / "S.bam", seed)
    write_stats(out / "S.stats", seed, insert_size="250")
    write_seed(d / "genome.fasta", seed)
    write_seed(d / "genome.gff", seed)

    r = run_arubaito("--steps", "genecounts", "rnaseq_pipeline", *other_settings(tmp_path),
                     f"out_dir={out}", f"bam_file={out}/S.bam", f"stats_file={out}/S.stats",
                     f"genome_file={d}/genome.fasta", f"gff_file={d}/genome.gff", env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    for path in [out / "S.genecounts.tab", out / "S.genecounts.tab.summary"]:
        check_output(path, seed, "featureCounts",
                     ["-T", "4", "-p", "--countReadPairs", "-s", "2", "-Q", "10", "-B", "-C", "-P",
                      "-d", "50", "-D", "1250", "-O", "-F", "GTF", "-t", "gene", "-g", "ID", "-f",
                      "--extraAttributes", "gene,gene_synonym,Dbxref,gene_biotype", "-a", f"{d}/genome.gff",
                      "--tmpDir", f"{tmp_path}/tmp/S-*", "-o", f"{out}/S.genecounts.tab", f"{out}/S.bam"],
                     inputs=[d / "genome.gff", out / "S.bam"])


def test_genecounts_nan_insert(tmp_path, seed, mock_env):
    d, out = tmp_path / "data", tmp_path / "data/P/results"
    write_seed(out / "S.bam", seed)
    write_stats(out / "S.stats", seed, insert_size="nan")
    write_seed(d / "genome.fasta", seed)
    write_seed(d / "genome.gff", seed)

    r = run_arubaito("--steps", "genecounts", "rnaseq_pipeline", *other_settings(tmp_path),
                     f"out_dir={out}", f"bam_file={out}/S.bam", f"stats_file={out}/S.stats",
                     f"genome_file={d}/genome.fasta", f"gff_file={d}/genome.gff", env=mock_env, cwd=tmp_path)
    assert r.returncode != 0

    assert "(got 'nan'), are the reads paired-end?" in read_text(out / "S.genecounts.log")
    assert not os.path.exists(out / "S.genecounts.tab")


def test_genecounts_missing_insert(tmp_path, seed, mock_env):
    d, out = tmp_path / "data", tmp_path / "data/P/results"
    write_seed(out / "S.bam", seed)
    write_seed(out / "S.stats", seed)
    write_seed(d / "genome.fasta", seed)
    write_seed(d / "genome.gff", seed)

    r = run_arubaito("--steps", "genecounts", "rnaseq_pipeline", *other_settings(tmp_path),
                     f"out_dir={out}", f"bam_file={out}/S.bam", f"stats_file={out}/S.stats",
                     f"genome_file={d}/genome.fasta", f"gff_file={d}/genome.gff", env=mock_env, cwd=tmp_path)
    assert r.returncode != 0

    assert "(got ''), are the reads paired-end?" in read_text(out / "S.genecounts.log")
    assert not os.path.exists(out / "S.genecounts.tab")


# ---------------------------------------------------------------------------
# refseq_pipeline

NCBI_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_000005845.2/download?"


def test_genome_download(tmp_path, seed, mock_env):
    refseq = tmp_path / "data/P/refseq"

    r = run_arubaito("--steps", "genome_download", "refseq_pipeline", *other_settings(tmp_path),
                     "refseq=GCF_000005845.2", env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    for path in [refseq / "GCF_000005845.2.fna", refseq / "GCF_000005845.2.gff"]:
        check_output(path, seed, "curl",
                     ["--fail", "-X", "GET",
                      NCBI_URL + "include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF"
                      "&hydrated=FULLY_HYDRATED",
                      "-H", "accept: application/zip", "-o", "ncbi_dataset.zip"])
    assert not os.path.exists(refseq / "GCF_000005845.2.gbff")


def test_genome_download_copy_outputs(tmp_path, seed, mock_env):
    refseq = tmp_path / "data/P/refseq"

    r = run_arubaito("--steps", "genome_download", "refseq_pipeline", *other_settings(tmp_path),
                     "refseq=GCF_000005845.2", "execute_copy_outputs=True", env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    for path in [refseq / "GCF_000005845.2.fna", refseq / "GCF_000005845.2.gff"]:
        check_output(path, seed, "curl",
                     ["--fail", "-X", "GET",
                      NCBI_URL + "include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF"
                      "&hydrated=FULLY_HYDRATED",
                      "-H", "accept: application/zip", "-o", "ncbi_dataset.zip"])


def test_genome_gbff_download(tmp_path, seed, mock_env):
    refseq = tmp_path / "data/P/refseq"

    r = run_arubaito("--steps", "genome_gbff_download", "refseq_pipeline", *other_settings(tmp_path),
                     "refseq=GCF_000005845.2", env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    check_output(refseq / "GCF_000005845.2.gbff", seed, "curl",
                 ["--fail", "-X", "GET",
                  NCBI_URL + "include_annotation_type=GENOME_GBFF&hydrated=FULLY_HYDRATED",
                  "-H", "accept: application/zip", "-o", "ncbi_dataset.zip"])
    assert not os.path.exists(refseq / "GCF_000005845.2.fna")


def test_genome_download_faonly(tmp_path, seed, mock_env):
    refseq = tmp_path / "data/P/refseq"

    r = run_arubaito("--steps", "genome_download_faonly", "refseq_pipeline", *other_settings(tmp_path),
                     "refseq=GCF_000005845.2", env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    check_output(refseq / "GCF_000005845.2.fna", seed, "curl",
                 ["--fail", "-X", "GET",
                  NCBI_URL + "include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF"
                  "&hydrated=FULLY_HYDRATED",
                  "-H", "accept: application/zip", "-o", "ncbi_dataset.zip"])
    assert not os.path.exists(refseq / "GCF_000005845.2.gff")


# ---------------------------------------------------------------------------
# sra_pipeline, sra_se_pipeline

def test_sra(tmp_path, seed, mock_env):
    sra = tmp_path / "data/P/sra"

    r = run_arubaito("--steps", "sra", "sra_pipeline", *other_settings(tmp_path), "accession=SRR123",
                     env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    check_output(sra / "SRR123.sra", seed, "prefetch",
                 ["-o", f"{sra}/SRR123.sra", "--resume", "no", "--verify", "yes", "--progress", "SRR123"])


def test_sra_single(tmp_path, seed, mock_env):
    sra = tmp_path / "data/P/sra"

    r = run_arubaito("--steps", "sra", "sra_se_pipeline", *other_settings(tmp_path), "accession=SRR123",
                     env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    check_output(sra / "SRR123.sra", seed, "prefetch",
                 ["-o", f"{sra}/SRR123.sra", "--resume", "no", "--verify", "yes", "--progress", "SRR123"])


def test_sra_fastq_paired_gz(tmp_path, seed, mock_env):
    sra, fq = tmp_path / "data/P/sra", tmp_path / "data/P/fastq"
    write_seed(sra / "SRR123.sra", seed)

    r = run_arubaito("--steps", "fastq", "sra_pipeline", *other_settings(tmp_path), "accession=SRR123",
                     f"read1_file={fq}/S_r1.fastq.gz", f"read2_file={fq}/S_r2.fastq.gz", env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    for path in [fq / "S_r1.fastq.gz", fq / "S_r2.fastq.gz"]:
        check_output(path, seed, "fastq-dump",
                     ["--outdir", f"{tmp_path}/tmp/S-*", "--split-3", "--log-level", "info", "-v", "-v",
                      f"{sra}/SRR123.sra"],
                     inputs=[sra / "SRR123.sra"])
    assert "Ignoring single-end reads" not in r.stdout


def test_sra_fastq_paired_uncompressed(tmp_path, seed, mock_env):
    sra, fq = tmp_path / "data/P/sra", tmp_path / "data/P/fastq"
    write_seed(sra / "SRR123.sra", seed)

    r = run_arubaito("--steps", "fastq", "sra_pipeline", *other_settings(tmp_path), "accession=SRR123",
                     f"read1_file={fq}/S_r1.fastq", f"read2_file={fq}/S_r2.fastq", env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    for path in [fq / "S_r1.fastq", fq / "S_r2.fastq"]:
        check_output(path, seed, "fastq-dump",
                     ["--outdir", f"{tmp_path}/tmp/S-*", "--split-3", "--log-level", "info", "-v", "-v",
                      f"{sra}/SRR123.sra"],
                     inputs=[sra / "SRR123.sra"])
        with open(path, "rb") as f:
            assert f.read(5) == b"tool\t"


def test_sra_fastq_paired_with_single_end_reads(tmp_path, seed, mock_env):
    sra, fq = tmp_path / "data/P/sra", tmp_path / "data/P/fastq"
    write_seed(sra / "SRR123.sra", seed)
    mock_env["MOCK_SRA_LAYOUT"] = "both"

    r = run_arubaito("--steps", "fastq", "sra_pipeline", *other_settings(tmp_path), "accession=SRR123",
                     f"read1_file={fq}/S_r1.fastq.gz", f"read2_file={fq}/S_r2.fastq.gz", env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    for path in [fq / "S_r1.fastq.gz", fq / "S_r2.fastq.gz"]:
        check_output(path, seed, "fastq-dump",
                     ["--outdir", f"{tmp_path}/tmp/S-*", "--split-3", "--log-level", "info", "-v", "-v",
                      f"{sra}/SRR123.sra"],
                     inputs=[sra / "SRR123.sra"])
    assert "WARNING: Ignoring single-end reads" in read_text(tmp_path / "data/P/results/S.fastq.log")


def test_sra_fastq_single_gz(tmp_path, seed, mock_env):
    sra, fq = tmp_path / "data/P/sra", tmp_path / "data/P/fastq"
    write_seed(sra / "SRR123.sra", seed)
    mock_env["MOCK_SRA_LAYOUT"] = "single"

    r = run_arubaito("--steps", "fastq", "sra_se_pipeline", *other_settings(tmp_path), "accession=SRR123",
                     f"read_file={fq}/S.fastq.gz", env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    check_output(fq / "S.fastq.gz", seed, "fastq-dump",
                 ["--outdir", f"{tmp_path}/tmp/S-*", "--split-3", "--log-level", "info", "-v", "-v",
                  f"{sra}/SRR123.sra"],
                 inputs=[sra / "SRR123.sra"])
    assert "Ignoring paired-end reads" not in r.stdout


def test_sra_fastq_single_uncompressed(tmp_path, seed, mock_env):
    sra, fq = tmp_path / "data/P/sra", tmp_path / "data/P/fastq"
    write_seed(sra / "SRR123.sra", seed)
    mock_env["MOCK_SRA_LAYOUT"] = "single"

    r = run_arubaito("--steps", "fastq", "sra_se_pipeline", *other_settings(tmp_path), "accession=SRR123",
                     f"read_file={fq}/S.fastq", env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    check_output(fq / "S.fastq", seed, "fastq-dump",
                 ["--outdir", f"{tmp_path}/tmp/S-*", "--split-3", "--log-level", "info", "-v", "-v",
                  f"{sra}/SRR123.sra"],
                 inputs=[sra / "SRR123.sra"])
    with open(fq / "S.fastq", "rb") as f:
        assert f.read(5) == b"tool\t"


def test_sra_fastq_single_with_paired_end_reads(tmp_path, seed, mock_env):
    sra, fq = tmp_path / "data/P/sra", tmp_path / "data/P/fastq"
    write_seed(sra / "SRR123.sra", seed)
    mock_env["MOCK_SRA_LAYOUT"] = "both"

    r = run_arubaito("--steps", "fastq", "sra_se_pipeline", *other_settings(tmp_path), "accession=SRR123",
                     f"read_file={fq}/S.fastq.gz", env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    check_output(fq / "S.fastq.gz", seed, "fastq-dump",
                 ["--outdir", f"{tmp_path}/tmp/S-*", "--split-3", "--log-level", "info", "-v", "-v",
                  f"{sra}/SRR123.sra"],
                 inputs=[sra / "SRR123.sra"])
    assert "WARNING: Ignoring paired-end reads" in read_text(tmp_path / "data/P/results/S.fastq.log")


# ---------------------------------------------------------------------------
# sra_pipeline2, sra_se_pipeline2

def test_sra2_fastq_paired_gz(tmp_path, seed, mock_env):
    fq = tmp_path / "data/P/fastq"

    r = run_arubaito("--steps", "fastq", "sra_pipeline2", *other_settings(tmp_path), "accession=SRR123",
                     f"read1_file={fq}/S_r1.fastq.gz", f"read2_file={fq}/S_r2.fastq.gz", env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    for path in [fq / "S_r1.fastq.gz", fq / "S_r2.fastq.gz"]:
        check_output(path, seed, "fasterq-dump", ["--split-3", "SRR123"])


def test_sra2_fastq_paired_uncompressed(tmp_path, seed, mock_env):
    fq = tmp_path / "data/P/fastq"

    r = run_arubaito("--steps", "fastq", "sra_pipeline2", *other_settings(tmp_path), "accession=SRR123",
                     f"read1_file={fq}/S_r1.fastq", f"read2_file={fq}/S_r2.fastq", env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    for path in [fq / "S_r1.fastq", fq / "S_r2.fastq"]:
        check_output(path, seed, "fasterq-dump", ["--split-3", "SRR123"])
        with open(path, "rb") as f:
            assert f.read(5) == b"tool\t"


def test_sra2_fastq_paired_several_accessions(tmp_path, seed, mock_env):
    fq = tmp_path / "data/P/fastq"

    r = run_arubaito("--steps", "fastq", "sra_pipeline2", *other_settings(tmp_path), "accession=SRR1,SRR2,SRR3",
                     f"read1_file={fq}/S_r1.fastq.gz", f"read2_file={fq}/S_r2.fastq.gz", env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    for path in [fq / "S_r1.fastq.gz", fq / "S_r2.fastq.gz"]:
        records = ["tool\t" + t for t in read_text(path).split("tool\t")[1:]]
        assert len(records) == 3
        check_record(records[0], seed, "fasterq-dump", ["--split-3", "SRR1"])
        check_record(records[1], seed, "fasterq-dump", ["--split-3", "SRR2"])
        check_record(records[2], seed, "fasterq-dump", ["--split-3", "SRR3"])


def test_sra2_fastq_paired_with_single_end_reads(tmp_path, seed, mock_env):
    fq = tmp_path / "data/P/fastq"
    mock_env["MOCK_SRA_LAYOUT"] = "both"

    r = run_arubaito("--steps", "fastq", "sra_pipeline2", *other_settings(tmp_path), "accession=SRR1,SRR2",
                     f"read1_file={fq}/S_r1.fastq.gz", f"read2_file={fq}/S_r2.fastq.gz", env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    log = read_text(tmp_path / "data/P/results/S.fastq.log")
    assert "WARNING: Ignoring single-end reads for accession SRR1" in log
    assert "WARNING: Ignoring single-end reads for accession SRR2" in log


def test_sra_se2_fastq_gz(tmp_path, seed, mock_env):
    fq = tmp_path / "data/P/fastq"
    mock_env["MOCK_SRA_LAYOUT"] = "single"

    r = run_arubaito("--steps", "fastq", "sra_se_pipeline2", *other_settings(tmp_path), "accession=SRR123",
                     f"read_file={fq}/S.fastq.gz", env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    check_output(fq / "S.fastq.gz", seed, "fasterq-dump", ["--split-3", "SRR123"])


def test_sra_se2_fastq_uncompressed(tmp_path, seed, mock_env):
    fq = tmp_path / "data/P/fastq"
    mock_env["MOCK_SRA_LAYOUT"] = "single"

    r = run_arubaito("--steps", "fastq", "sra_se_pipeline2", *other_settings(tmp_path), "accession=SRR123",
                     f"read_file={fq}/S.fastq", env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    check_output(fq / "S.fastq", seed, "fasterq-dump", ["--split-3", "SRR123"])
    with open(fq / "S.fastq", "rb") as f:
        assert f.read(5) == b"tool\t"


def test_sra_se2_fastq_several_accessions(tmp_path, seed, mock_env):
    fq = tmp_path / "data/P/fastq"
    mock_env["MOCK_SRA_LAYOUT"] = "single"

    r = run_arubaito("--steps", "fastq", "sra_se_pipeline2", *other_settings(tmp_path), "accession=SRR1,SRR2",
                     f"read_file={fq}/S.fastq.gz", env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    records = ["tool\t" + t for t in read_text(fq / "S.fastq.gz").split("tool\t")[1:]]
    assert len(records) == 2
    check_record(records[0], seed, "fasterq-dump", ["--split-3", "SRR1"])
    check_record(records[1], seed, "fasterq-dump", ["--split-3", "SRR2"])


def test_sra_se2_fastq_with_paired_end_reads(tmp_path, seed, mock_env):
    fq = tmp_path / "data/P/fastq"
    mock_env["MOCK_SRA_LAYOUT"] = "both"

    r = run_arubaito("--steps", "fastq", "sra_se_pipeline2", *other_settings(tmp_path), "accession=SRR1",
                     f"read_file={fq}/S.fastq.gz", env=mock_env, cwd=tmp_path)
    assert r.returncode == 0

    assert "WARNING: Ignoring paired-end reads for accession SRR1" in read_text(
        tmp_path / "data/P/results/S.fastq.log")
