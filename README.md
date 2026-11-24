# MFASeqPipeline

## Overview

`MFASeqPipeline` is a simple but flexible pipeline for analysing *marker
frequency analysis* (MFA-seq) data. Its first stage maps sequencing reads and
produces a table of position-wise DNA abundances, comparing replicating
(exponentially growing) samples with non-replicating (e.g. stationary)
references. From these abundances, it then infers replication dynamics, in
particular the speed of replication forks over time with
[VeloGP](https://github.com/oist/VeloGP). The main features are

  - Non-parametric, Gaussian-process based estimation of fork velocity
  - Automatic detection and exclusion of genomic regions that cannot be mapped
    reliably
  - Works for arbitrary organisms, as long as a reference genome is available
  - Runs locally or on an HPC cluster (Slurm by default)
  - Companion pipelines to download data from SRA and genomes from NCBI, and to
    count RNA-seq reads per gene

The pipelines are built on **arubaito**, a simple workflow engine contained in
the script `arubaito.py` and installed as the command `arubaito`. Pipelines are
YAML files that define processing steps as shell scripts, and a project is a
small YAML file that names the samples and the steps to run for each. See [docs/arubaito.md](docs/arubaito.md).

## Installation

MFASeqPipeline is installed with pip into a conda environment, which provides
the external tools the pipelines use (bwa, samtools, sra-tools, Subread, R,
PyTorch, gnuplot, …). First download the environment file and create and
activate the environment:

```bash
curl -O https://raw.githubusercontent.com/oist/MFASeqPipeline/master/env.yaml
```

```bash
conda env create -f env.yaml
```

```bash
conda activate mfa_pipeline
```

Then install MFASeqPipeline into it:

```bash
pip install git+https://github.com/oist/MFASeqPipeline.git
```

This installs the command `arubaito` together with all pipelines. The pipelines
are found inside the installed package, so project files can live in any
directory. Run `arubaito` from the directory that contains your project file.

### Alternative: running from a clone without installing

Instead of installing with pip, you can also check out the repository, create
the environment from it, and run `arubaito.py` from the checkout directly. Clone
the repository together with its submodule VeloGP:

```bash
git clone --recurse-submodules git@github.com:oist/MFASeqPipeline.git
```

```bash
conda env create -f MFASeqPipeline/env.yaml
```

```bash
conda activate mfa_pipeline
```

Then use `/path/to/MFASeqPipeline/arubaito.py` wherever this documentation uses
the `arubaito` command, for example:

```bash
/path/to/MFASeqPipeline/arubaito.py myproject genome stationary exponential
```

If you download the repository as an archive instead of cloning it, the archive
does not contain VeloGP; download [VeloGP](https://github.com/oist/VeloGP) as
well and put it into `MFASeqPipeline/steps/velogp`.

To get the `arubaito` command while still using the files in the checkout (e.g.
to modify the pipelines), install the checkout in editable mode:

```bash
pip install -e MFASeqPipeline
```

## Quick start

This example analyses two stationary-phase (non-replicating) samples, `stat1`
and `stat2`, and two exponentially growing samples, `exp1` and `exp2`. Each
growing sample is compared with one stationary sample. The reads are
paired-end, and everything runs on the local machine.

Arrange the input files like this:

```
/data/genome.fasta                       # reference genome (FASTA)
/data/genome.gff                         # genome annotation (GFF3)
/data/origins.bed                        # replication origin(s), named origin_of_replication
/data/myproject/fastq/stat1_r1.fastq.gz  /data/myproject/fastq/stat1_r2.fastq.gz
/data/myproject/fastq/stat2_r1.fastq.gz  /data/myproject/fastq/stat2_r2.fastq.gz
/data/myproject/fastq/exp1_r1.fastq.gz   /data/myproject/fastq/exp1_r2.fastq.gz
/data/myproject/fastq/exp2_r1.fastq.gz   /data/myproject/fastq/exp2_r2.fastq.gz
/data/myproject/meta/growthrates.tab     # growth rates of exp1 and exp2
```

`growthrates.tab` lists each growing sample's growth rate in 1/hour, tab-separated
(see [docs/growthrates_format.md](docs/growthrates_format.md)):

```
exp1	1.386
exp2	0.693
```

Create the project file `myproject.yaml`:

```yaml
use: mfa_pipeline
project: myproject
base_dir: /data                 # project directory is {base_dir}/{project}
# Run steps on this machine instead of submitting them with Slurm
execute_cmd: '{script_file}'
execute_copy_outputs: False
cores: 4                        # cores per step
targets:
  genome:
    - target_label: genome
      steps: [ genome_index ]
  stationary:
    override:
      steps:
        - map
        - fastqc
        - mapstats
        - mappability
    targets:
      - sample: stat1
      - sample: stat2
  exponential:
    override:
      steps:
        - map
        - fastqc
        - mapstats
        - mappability
        - abundance
        - velogp
    targets:
      - sample: exp1
        reference: stat1
      - sample: exp2
        reference: stat2
```

Run the three targets in order. The genome index comes first, and the stationary
samples come before the exponential samples that use them as references. `--jobs 2`
processes two samples at a time:

```bash
arubaito --jobs 2 myproject genome stationary exponential
```

Results are written to `/data/myproject/results/`. The main ones are:

| File | Content |
|---|---|
| `<sample>.bam` | Mapped reads |
| `<sample>.mappability.bed` | Regions that cannot be mapped reliably |
| `exp1.abd.w10kbp.gz` | Abundances in 10 kb windows, relative to `stat1` ([format](docs/abundance_format.md)) |
| `exp1.velogp.vt.gz`, `.vt.pdf` | Fork velocity over time, estimated by VeloGP |
| `<sample>.<step>.log` | Log of each step |

Steps whose outputs already exist are skipped, so you can rerun the same command
after adding samples or fixing a failed step. To run on a Slurm cluster, remove
the `execute_cmd` and `execute_copy_outputs` lines (see
[HPC clusters](docs/arubaito.md#hpc-clusters)).

## Documentation

* [docs/arubaito.md](docs/arubaito.md): how arubaito works: steps, variables,
  targets, and running on HPC clusters.
* [docs/summary.md](docs/summary.md): a compact technical reference of the engine.

### Pipelines

| Pipeline | Purpose |
|---|---|
| [mfa_pipeline](docs/mfa_pipeline.md) | MFA-seq analysis of paired-end data: mapping, mappability, abundances, VeloGP and VeloMix |
| [mfa_se_pipeline](docs/mfa_se_pipeline.md) | The same for single-end data |
| [sra_pipeline](docs/sra_pipeline.md), [sra_se_pipeline](docs/sra_se_pipeline.md) | Download paired-end / single-end runs from SRA (via `.sra` files), for use with `mfa_pipeline` / `mfa_se_pipeline` |
| [sra_pipeline2](docs/sra_pipeline2.md), [sra_se_pipeline2](docs/sra_se_pipeline2.md) | Download paired-end / single-end runs from SRA directly to FASTQ, several runs per sample |
| [refseq_pipeline](docs/refseq_pipeline.md) | Download a reference genome and annotation from NCBI |
| [rnaseq_pipeline](docs/rnaseq_pipeline.md) | Count RNA-seq fragments per gene, for use with `mfa_pipeline` |

The overlay `mfa_no_reference.yaml` adapts `mfa_pipeline` to samples without a
reference; see [abundance](docs/mfa_pipeline.md#abundance).

### File formats

* [Abundance file](docs/abundance_format.md): output of `abundance`
* [Mappability file](docs/mappability_format.md): output of `mappability`
* [Mapping statistics](docs/stats_format.md): output of `mapstats`
* [Growth rates](docs/growthrates_format.md): input for `velogp`
* [Velocity oscillation fit](docs/velocity_oscillation_format.md): output of the
  outdated `velocity_oscillation` step
