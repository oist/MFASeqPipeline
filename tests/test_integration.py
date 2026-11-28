"""Tests that run whole pipelines with arubaito, using the mock tools.

Each test writes a project file and the workflow's input files (containing just
the seed), runs the targets, and checks every output: that it was written by the
expected tool, with the expected arguments, from the expected inputs. Inputs are
the files produced by earlier steps, so the checks also verify that the steps are
chained correctly.

Test names must not contain "insert_size" or "read_length", see test_steps.py.
"""
import os

import yaml

from conftest import (MOCK_PIPELINE, check_output, check_record, read_text, run_arubaito,
                      write_growthrates, write_seed)

NCBI_URL = ("https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_1/download?"
            "include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&hydrated=FULLY_HYDRATED")


def test_paired_end_from_sra_and_refseq(tmp_path, seed, mock_env):
    """Download a genome and two paired-end samples, then run all MFA steps. Inputs
    and outputs are copied to and from the temporary directory, as on a cluster."""
    d, out = tmp_path / "data", tmp_path / "data/P/results"
    refseq, sra, fq = tmp_path / "data/P/refseq", tmp_path / "data/P/sra", tmp_path / "data/P/fastq"
    G, T, Ts = f"{tmp_path}/tmp/P/genome-*", f"{tmp_path}/tmp/P/exp-*", f"{tmp_path}/tmp/P/stat-*"
    write_growthrates(d / "P/meta/growthrates.tab", seed, "exp", "0.693")
    (tmp_path / "testproject.yaml").write_text(yaml.safe_dump({
        "use": ["sra_pipeline", "refseq_pipeline", "mfa_pipeline"],
        "mfa_pipeline_dir": MOCK_PIPELINE,
        "base_dir": str(d),
        "project": "P",
        "tmp_base_dir": f"{tmp_path}/tmp",
        "execute_cmd": "{script_file}",
        "execute_copy_inputs": True,
        "execute_copy_outputs": True,
        "cores": 2,
        "refseq": "GCF_1",
        "targets": {
            "genome": [{"target_label": "genome", "steps": ["genome_download", "genome_index"]}],
            "stationary": [{"sample": "stat", "accession": "SRR1",
                            "steps": ["sra", "fastq", "map", "fastqc", "mapstats", "mappability"]}],
            "exponential": [{"sample": "exp", "accession": "SRR2", "reference": "stat",
                             "steps": ["sra", "fastq", "map", "fastqc", "mapstats", "mappability", "abundance",
                                       "abundance_plot", "velogp", "velomix"]}],
        },
    }))

    r = run_arubaito("testproject", "genome", "stationary", "exponential", env=mock_env, cwd=tmp_path)
    assert r.returncode == 0
    assert "JOB FAILED" not in r.stdout

    # genome
    for path in [refseq / "GCF_1.fna", refseq / "GCF_1.gff"]:
        check_output(path, seed, "curl",
                     ["--fail", "-X", "GET", NCBI_URL, "-H", "accept: application/zip", "-o", "ncbi_dataset.zip"])
    check_output(refseq / "GCF_1.fna.bwaidx/index.mockidx", seed, "bwa",
                 ["index", "-p", f"{G}/*/GCF_1.fna.bwaidx/index", f"{G}/*/GCF_1.fna"],
                 inputs=[refseq / "GCF_1.fna"])

    # stationary and exponential samples
    for sample, accession, tmp in [("stat", "SRR1", Ts), ("exp", "SRR2", T)]:
        check_output(sra / f"{accession}.sra", seed, "prefetch",
                     ["-o", f"{tmp}/*/{accession}.sra", "--resume", "no", "--verify", "yes", "--progress", accession])
        for path in [fq / f"{sample}_r1.fastq.gz", fq / f"{sample}_r2.fastq.gz"]:
            check_output(path, seed, "fastq-dump",
                         ["--outdir", tmp, "--split-3", "--log-level", "info", "-v", "-v", f"{tmp}/*/{accession}.sra"],
                         inputs=[sra / f"{accession}.sra"])
        check_output(out / f"{sample}.bam", seed, "samtools",
                     ["sort", "--threads", "2", "-m", "8192M", "-T", f"{tmp}/{sample}.samtools-sort.", "-O", "BAM",
                      "-l", "7", "-o", f"{tmp}/*/{sample}.bam.in-progress", "-"],
                     stdin=dict(tool="bwa",
                                args=["mem", "-t", "2", "-v", "1", f"{tmp}/*/GCF_1.fna.bwaidx/index",
                                      f"{tmp}/{sample}.r1.pipe", f"{tmp}/{sample}.r2.pipe"],
                                inputs=[refseq / "GCF_1.fna.bwaidx/index.mockidx",
                                        (f"{tmp}/{sample}.r1.pipe", fq / f"{sample}_r1.fastq.gz"),
                                        (f"{tmp}/{sample}.r2.pipe", fq / f"{sample}_r2.fastq.gz")]))
        check_output(out / f"{sample}.bai", seed, "samtools",
                     ["index", "-@", "2", f"{tmp}/*/{sample}.bam", f"{tmp}/*/{sample}.bai"],
                     inputs=[out / f"{sample}.bam"])
        for path in [out / f"{sample}_fastqc.html", out / f"{sample}_fastqc.zip"]:
            check_output(path, seed, "fastqc", ["--threads", "2", "--noextract", f"{tmp}/*/{sample}.bam"],
                         inputs=[out / f"{sample}.bam"])
        check_output(out / f"{sample}.stats", seed, "mapping_statistics.py",
                     [f"{tmp}/*/{sample}.bam", f"{tmp}/*/{sample}.stats"],
                     inputs=[out / f"{sample}.bam"])
        check_output(out / f"{sample}.mappability.bed", seed, "mappability.py",
                     ["--genome", f"{tmp}/*/GCF_1.fna", "--genome-index", f"{tmp}/*/GCF_1.fna.bwaidx/index",
                      "--circular", "--paired", "--fragment-length", "250", "--read-length", "100",
                      "--mapping-quality", "20", "--regions", "", "--exclude-regions", "",
                      "--output", f"{tmp}/*/{sample}.mappability.bed", "--", "-t", "2", "-v", "1"],
                     inputs=[refseq / "GCF_1.fna", refseq / "GCF_1.fna.bwaidx/index.mockidx"])

    # exponential sample only
    check_output(out / "exp.abd.w10kbp.gz", seed, "abundance.py",
                 ["--sample-bam", f"{T}/*/exp.bam", "--sample-mappability", f"{T}/*/exp.mappability.bed",
                  "--reference-bam", f"{T}/*/stat.bam", "--reference-mappability", f"{T}/*/stat.mappability.bed",
                  "--paired", "--mapping-quality", "20", "--windowsize", "10k",
                  "--output", f"{T}/*/exp.abd.w10kbp.gz"],
                 inputs=[out / "exp.bam", out / "exp.mappability.bed", out / "stat.bam", out / "stat.mappability.bed"])
    check_output(out / "exp.abd.w10kbp.pdf", seed, "gnuplot", [],
                 inputs=[out / "exp.abd.w10kbp.gz"], stdin_contains=['column("coverage_ratio")'])
    for name in ["abd.gz", "abd.pdf", "vt.gz", "vt.pdf", "vx.gz", "pkl"]:
        check_output(out / f"exp.velogp.{name}", seed, "velogp",
                     ["--cores", "4", "--genome", f"{T}/*/GCF_1.fna", "--gff", f"{T}/*/GCF_1.gff",
                      "--columns", "seq,midpoint,coverage_ratio,coverage_ratio_stderr,ratio_mappable",
                      "--vx-feature-types", "", "--vx-feature-names", "", "--vx-delayscale", "0.01",
                      "--vt-covariance", "m1", "--vt-timescale", "0.2", "--vt-magnitude", "0.25",
                      "--out-abundance", f"{T}/*/exp.velogp.abd.gz", "--out-trajectories", f"{T}/*/exp.velogp.vt.gz",
                      "--out-vx-feature-delays", f"{T}/*/exp.velogp.vx.gz",
                      "--out-objects", f"{T}/*/exp.velogp.pkl",
                      "--plot-vt", f"{T}/*/exp.velogp.vt.pdf", "--plot-abundance", f"{T}/*/exp.velogp.abd.pdf",
                      "--origins", "origin_of_replication", "--growthrate", "0.693", "--seed", "42",
                      "--perturbation-rejects", "5", "--perturbation-reduction", "0.8408964",
                      "--perturbation-min", "0.01", f"{T}/*/exp.abd.w10kbp.gz"],
                     inputs=[out / "exp.abd.w10kbp.gz", refseq / "GCF_1.fna", refseq / "GCF_1.gff"])
    for name in ["abd.gz", "abd.pdf", "p.gz", "vt.gz", "vt.pdf", "vx.gz", "pkl"]:
        check_output(out / f"exp.velomix.{name}", seed, "velogp.velomix",
                     ["--stall-sites", "", "--cores", "4",
                      "--out-probabilities", f"{T}/*/exp.velomix.p.gz",
                      "--out-abundance", f"{T}/*/exp.velomix.abd.gz",
                      "--out-trajectories", f"{T}/*/exp.velomix.vt.gz",
                      "--out-vx-feature-delays", f"{T}/*/exp.velomix.vx.gz",
                      "--out-objects", f"{T}/*/exp.velomix.pkl",
                      "--plot-vt", f"{T}/*/exp.velomix.vt.pdf", "--plot-abundance", f"{T}/*/exp.velomix.abd.pdf",
                      "--seed", "42", "--perturbation-rejects", "5", "--perturbation-reduction", "0.8408964",
                      "--perturbation-min", "0.01", f"{T}/*/exp.velogp.pkl"],
                     inputs=[out / "exp.velogp.pkl"])

    # Logs, and no leftover temporary directories
    for log in ["GCF_1.genome_download", "genome.genome_index", "exp.sra", "exp.fastq", "exp.map", "exp.fastqc",
                "exp.mapstats", "exp.mappability", "exp.abundance.w10kbp", "exp.abundance_plot", "exp.velogp",
                "exp.velomix", "stat.map"]:
        assert os.path.isfile(out / f"{log}.log"), f"log {log}.log missing"
    assert os.listdir(tmp_path / "tmp/P") == []


def test_single_end_without_reference(tmp_path, seed, mock_env):
    """Download single-end samples from several runs each and analyse them without a reference,
    running steps locally (no copying)."""
    d, out, fq = tmp_path / "data", tmp_path / "data/P/results", tmp_path / "data/P/fastq"
    T = f"{tmp_path}/tmp/P/S-*"
    write_seed(d / "genome.fasta", seed)
    write_seed(d / "genome.gff", seed)
    write_growthrates(d / "P/meta/growthrates.tab", seed, "S", "1.5")
    mock_env["MOCK_SRA_LAYOUT"] = "single"
    (tmp_path / "testproject.yaml").write_text(yaml.safe_dump({
        "use": ["sra_se_pipeline2", "mfa_no_reference", "mfa_se_pipeline"],
        "mfa_pipeline_dir": MOCK_PIPELINE,
        "base_dir": str(d),
        "project": "P",
        "tmp_base_dir": f"{tmp_path}/tmp",
        "execute_cmd": "{script_file}",
        "execute_copy_outputs": False,
        "cores": 4,
        "genome_circular": False,
        "windowsize": "5k",
        "targets": {
            "genome": [{"target_label": "genome", "steps": ["genome_index"]}],
            "samples": [{"sample": "S", "accession": "SRR1,SRR2",
                         "steps": ["fastq", "map", "mapstats", "mappability", "abundance", "abundance_plot",
                                   "velogp"]}],
        },
    }))

    r = run_arubaito("testproject", "genome", "samples", env=mock_env, cwd=tmp_path)
    assert r.returncode == 0
    assert "JOB FAILED" not in r.stdout

    check_output(d / "genome.fasta.bwaidx/index.mockidx", seed, "bwa",
                 ["index", "-p", f"{d}/genome.fasta.bwaidx/index", f"{d}/genome.fasta"],
                 inputs=[d / "genome.fasta"])
    records = ["tool\t" + t for t in read_text(fq / "S.fastq.gz").split("tool\t")[1:]]
    assert len(records) == 2
    check_record(records[0], seed, "fasterq-dump", ["--split-3", "SRR1"])
    check_record(records[1], seed, "fasterq-dump", ["--split-3", "SRR2"])
    check_output(out / "S.bam", seed, "samtools",
                 ["sort", "--threads", "4", "-m", "4096M", "-T", f"{T}/S.samtools-sort.", "-O", "BAM",
                  "-l", "7", "-o", f"{out}/S.bam.in-progress", "-"],
                 stdin=dict(tool="bwa",
                            args=["mem", "-t", "4", "-v", "1", f"{d}/genome.fasta.bwaidx/index", f"{T}/S.r.pipe"],
                            inputs=[d / "genome.fasta.bwaidx/index.mockidx", (f"{T}/S.r.pipe", fq / "S.fastq.gz")]))
    check_output(out / "S.bai", seed, "samtools", ["index", "-@", "4", f"{out}/S.bam", f"{out}/S.bai"],
                 inputs=[out / "S.bam"])
    check_output(out / "S.stats", seed, "mapping_statistics.py", [f"{out}/S.bam", f"{out}/S.stats"],
                 inputs=[out / "S.bam"])
    check_output(out / "S.mappability.bed", seed, "mappability.py",
                 ["--genome", f"{d}/genome.fasta", "--genome-index", f"{d}/genome.fasta.bwaidx/index",
                  "--read-length", "100", "--mapping-quality", "20", "--regions", "", "--exclude-regions", "",
                  "--output", f"{out}/S.mappability.bed", "--", "-t", "4", "-v", "1"],
                 inputs=[d / "genome.fasta", d / "genome.fasta.bwaidx/index.mockidx"])
    check_output(out / "S.abd.w5kbp.gz", seed, "abundance.py",
                 ["--sample-bam", f"{out}/S.bam", "--sample-mappability", f"{out}/S.mappability.bed",
                  "--mapping-quality", "20", "--windowsize", "5k", "--output", f"{out}/S.abd.w5kbp.gz"],
                 inputs=[out / "S.bam", out / "S.mappability.bed"])
    check_output(out / "S.abd.w5kbp.pdf", seed, "gnuplot", [],
                 inputs=[out / "S.abd.w5kbp.gz"], stdin_contains=['column("sample_coverage_avg")'])
    for name in ["abd.gz", "abd.pdf", "vt.gz", "vt.pdf", "vx.gz", "pkl"]:
        check_output(out / f"S.velogp.{name}", seed, "velogp",
                     ["--cores", "4", "--genome", f"{d}/genome.fasta", "--gff", f"{d}/genome.gff",
                      "--columns", "seq,midpoint,sample_coverage_avg,sample_coverage_avg_stderr,sample_mappable",
                      "--vx-feature-types", "", "--vx-feature-names", "", "--vx-delayscale", "0.01",
                      "--vt-covariance", "m1", "--vt-timescale", "0.2", "--vt-magnitude", "0.25",
                      "--out-abundance", f"{out}/S.velogp.abd.gz", "--out-trajectories", f"{out}/S.velogp.vt.gz",
                      "--out-vx-feature-delays", f"{out}/S.velogp.vx.gz", "--out-objects", f"{out}/S.velogp.pkl",
                      "--plot-vt", f"{out}/S.velogp.vt.pdf", "--plot-abundance", f"{out}/S.velogp.abd.pdf",
                      "--origins", "origin_of_replication", "--growthrate", "1.5", "--seed", "42",
                      "--perturbation-rejects", "5", "--perturbation-reduction", "0.8408964",
                      "--perturbation-min", "0.01", f"{out}/S.abd.w5kbp.gz"],
                     inputs=[out / "S.abd.w5kbp.gz", d / "genome.fasta", d / "genome.gff"])


def test_rnaseq(tmp_path, seed, mock_env):
    """Map an RNA-seq sample with mfa_pipeline and count fragments per gene."""
    d, out, fq = tmp_path / "data", tmp_path / "data/P/results", tmp_path / "data/P/fastq"
    T = f"{tmp_path}/tmp/P/rna-*"
    write_seed(d / "genome.fasta", seed)
    write_seed(d / "genome.gff", seed)
    write_seed(fq / "rna_r1.fastq.gz", seed)
    write_seed(fq / "rna_r2.fastq.gz", seed)
    (tmp_path / "testproject.yaml").write_text(yaml.safe_dump({
        "use": ["rnaseq_pipeline", "mfa_pipeline"],
        "mfa_pipeline_dir": MOCK_PIPELINE,
        "base_dir": str(d),
        "project": "P",
        "tmp_base_dir": f"{tmp_path}/tmp",
        "execute_cmd": "{script_file}",
        "execute_copy_outputs": False,
        "cores": 4,
        "targets": {
            "genome": [{"target_label": "genome", "steps": ["genome_index"]}],
            "rnaseq": [{"sample": "rna", "steps": ["map", "mapstats", "genecounts"]}],
        },
    }))

    r = run_arubaito("testproject", "genome", "rnaseq", env=mock_env, cwd=tmp_path)
    assert r.returncode == 0
    assert "JOB FAILED" not in r.stdout

    check_output(out / "rna.bam", seed, "samtools",
                 ["sort", "--threads", "4", "-m", "4096M", "-T", f"{T}/rna.samtools-sort.", "-O", "BAM",
                  "-l", "7", "-o", f"{out}/rna.bam.in-progress", "-"],
                 stdin=dict(tool="bwa",
                            args=["mem", "-t", "4", "-v", "1", f"{d}/genome.fasta.bwaidx/index",
                                  f"{T}/rna.r1.pipe", f"{T}/rna.r2.pipe"],
                            inputs=[d / "genome.fasta.bwaidx/index.mockidx",
                                    (f"{T}/rna.r1.pipe", fq / "rna_r1.fastq.gz"),
                                    (f"{T}/rna.r2.pipe", fq / "rna_r2.fastq.gz")]))
    check_output(out / "rna.stats", seed, "mapping_statistics.py", [f"{out}/rna.bam", f"{out}/rna.stats"],
                 inputs=[out / "rna.bam"])
    for path in [out / "rna.genecounts.tab", out / "rna.genecounts.tab.summary"]:
        check_output(path, seed, "featureCounts",
                     ["-T", "4", "-p", "--countReadPairs", "-s", "2", "-Q", "10", "-B", "-C", "-P",
                      "-d", "50", "-D", "1250", "-O", "-F", "GTF", "-t", "gene", "-g", "ID", "-f",
                      "--extraAttributes", "gene,gene_synonym,Dbxref,gene_biotype", "-a", f"{d}/genome.gff",
                      "--tmpDir", T, "-o", f"{out}/rna.genecounts.tab", f"{out}/rna.bam"],
                     inputs=[d / "genome.gff", out / "rna.bam"])


def test_example_dataset(tmp_path, seed, mock_env):
    """Run the shipped example_dataset.yaml, with its placeholder paths replaced on the command line."""
    data, out, meta = tmp_path / "data", tmp_path / "results", tmp_path / "meta"
    write_seed(data / "genome.fasta", seed)
    write_seed(data / "genome.gff", seed)
    write_growthrates(meta / "growthrates.tab", seed, "exponential_sample_1", "0.7")
    with open(meta / "growthrates.tab", "a") as f:
        f.write("exponential_sample_2\t1.4\n")
    samples = ["stationary_sample_1", "stationary_sample_2", "exponential_sample_1", "exponential_sample_2"]
    for sample in samples:
        write_seed(data / f"{sample}.fastq.gz", seed)

    r = run_arubaito("example_dataset", "genome", "stationary", "exponential",
                     f"mfa_pipeline_dir={MOCK_PIPELINE}", f"read_file='{data}/{{sample}}.fastq.gz'",
                     f"meta_dir={meta}", f"out_dir={out}", f"tmp_dir='{tmp_path}/tmp/{{target_label}}-'",
                     f"genome_file={data}/genome.fasta", f"gff_file={data}/genome.gff",
                     env=mock_env, cwd=tmp_path)
    assert r.returncode == 0
    assert "JOB FAILED" not in r.stdout

    check_output(data / "genome.fasta.bwaidx/index.mockidx", seed, "bwa",
                 ["index", "-p", f"{data}/genome.fasta.bwaidx/index", f"{data}/genome.fasta"],
                 inputs=[data / "genome.fasta"])
    for sample in samples:
        T = f"{tmp_path}/tmp/{sample}-*"
        check_output(out / f"{sample}.bam", seed, "samtools",
                     ["sort", "--threads", "32", "-m", "512M", "-T", f"{T}/{sample}.samtools-sort.", "-O", "BAM",
                      "-l", "7", "-o", f"{out}/{sample}.bam.in-progress", "-"],
                     stdin=dict(tool="bwa",
                                args=["mem", "-t", "32", "-v", "1", f"{data}/genome.fasta.bwaidx/index",
                                      f"{T}/{sample}.r.pipe"],
                                inputs=[data / "genome.fasta.bwaidx/index.mockidx",
                                        (f"{T}/{sample}.r.pipe", data / f"{sample}.fastq.gz")]))
        check_output(out / f"{sample}.stats", seed, "mapping_statistics.py",
                     [f"{out}/{sample}.bam", f"{out}/{sample}.stats"], inputs=[out / f"{sample}.bam"])
        check_output(out / f"{sample}.mappability.bed", seed, "mappability.py",
                     ["--genome", f"{data}/genome.fasta", "--genome-index", f"{data}/genome.fasta.bwaidx/index",
                      "--circular", "--read-length", "100", "--mapping-quality", "20",
                      "--regions", "", "--exclude-regions", "",
                      "--output", f"{out}/{sample}.mappability.bed", "--", "-t", "32", "-v", "1"],
                     inputs=[data / "genome.fasta", data / "genome.fasta.bwaidx/index.mockidx"])
    for sample, reference, rate in [("exponential_sample_1", "stationary_sample_1", "0.7"),
                                    ("exponential_sample_2", "stationary_sample_2", "1.4")]:
        check_output(out / f"{sample}.abd.w10kbp.gz", seed, "abundance.py",
                     ["--sample-bam", f"{out}/{sample}.bam", "--sample-mappability", f"{out}/{sample}.mappability.bed",
                      "--reference-bam", f"{out}/{reference}.bam",
                      "--reference-mappability", f"{out}/{reference}.mappability.bed",
                      "--mapping-quality", "20", "--windowsize", "10k", "--output", f"{out}/{sample}.abd.w10kbp.gz"],
                     inputs=[out / f"{sample}.bam", out / f"{sample}.mappability.bed",
                             out / f"{reference}.bam", out / f"{reference}.mappability.bed"])
        rec = check_output(out / f"{sample}.velogp.pkl", seed, "velogp",
                           ["--cores", "4", "--genome", f"{data}/genome.fasta", "--gff", f"{data}/genome.gff",
                            "--columns", "seq,midpoint,coverage_ratio,coverage_ratio_stderr,ratio_mappable",
                            "--vx-feature-types", "", "--vx-feature-names", "", "--vx-delayscale", "0.01",
                            "--vt-covariance", "m1", "--vt-timescale", "0.2", "--vt-magnitude", "0.25",
                            "--out-abundance", f"{out}/{sample}.velogp.abd.gz",
                            "--out-trajectories", f"{out}/{sample}.velogp.vt.gz",
                            "--out-vx-feature-delays", f"{out}/{sample}.velogp.vx.gz",
                            "--out-objects", f"{out}/{sample}.velogp.pkl",
                            "--plot-vt", f"{out}/{sample}.velogp.vt.pdf",
                            "--plot-abundance", f"{out}/{sample}.velogp.abd.pdf",
                            "--origins", "origin_of_replication", "--growthrate", rate, "--seed", "42",
                            "--perturbation-rejects", "5", "--perturbation-reduction", "0.8408964",
                            "--perturbation-min", "0.01", f"{out}/{sample}.abd.w10kbp.gz"],
                           inputs=[out / f"{sample}.abd.w10kbp.gz", data / "genome.fasta", data / "genome.gff"])
        for name in ["abd.gz", "abd.pdf", "vt.gz", "vt.pdf", "vx.gz"]:
            assert read_text(out / f"{sample}.velogp.{name}") == read_text(out / f"{sample}.velogp.pkl")
