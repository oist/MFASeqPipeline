# mfa_pipeline

File: `mfa_pipeline.yaml`. Single-end variant:
[mfa_se_pipeline](mfa_se_pipeline.md). For how configurations, steps and
targets work in general, see [summary.md](summary.md).

## Summary

`mfa_pipeline` analyses **marker frequency analysis (MFA-seq)** data from
paired-end sequencing. It maps reads to a reference genome and works out which
genomic regions can be mapped reliably. It then computes position-wise DNA
abundances in windows, usually as the ratio between a replicating
(exponentially growing) sample and a non-replicating reference (e.g. a
stationary-phase sample). From this abundance profile, the downstream steps
infer replication dynamics using Gaussian-process models of fork speed
(`velogp`, `velomix`). The older parametric fit
`velocity_oscillation` is still available but outdated.

For samples without a reference, use the overlay `mfa_no_reference.yaml`; see
[abundance](#abundance).

## Quick start

Assume a directory `/data/myproject` holding a reference genome and gzipped
paired-end FASTQ files:

```
/data/genome.fasta
/data/genome.gff
/data/myproject/fastq/stat_r1.fastq.gz   /data/myproject/fastq/stat_r2.fastq.gz
/data/myproject/fastq/exp_r1.fastq.gz    /data/myproject/fastq/exp_r2.fastq.gz
/data/myproject/meta/growthrates.tab     # see growthrates_format.md
```

Create `myproject.yaml`:

```yaml
use: mfa_pipeline
project: myproject
base_dir: /data                 # project_dir = {base_dir}/{project}
# Run steps locally instead of with Slurm
execute_cmd: '{script_file}'
execute_copy_outputs: False
cores: 4
targets:
  genome:
    - target_label: genome
      steps: [ genome_index ]
  reference:
    - sample: stat
      steps: [ map, mapstats, mappability ]
  growing:
    - sample: exp
      reference: stat
      steps: [ map, fastqc, mapstats, mappability, abundance, abundance_plot, velogp ]
```

velogp looks up the replication origins in `genome.gff`, by default all
features of type `origin_of_replication`. Set `origins` to select them by name
or feature type instead, e.g. `origins: oriC`. If the GFF file does not contain
the origins, list them in a BED file and set `origin_file` to it (see [velogp](#velogp)).

Run the targets in dependency order:

```bash
arubaito myproject genome reference growing
```

Results go to `/data/myproject/results/`. For example,
`exp.abd.w10kbp.gz` holds the abundances and `exp.velogp.vt.gz` the fitted fork
speed over time. Each step writes a log file, `<sample>.<step>.log`, next to its results.
To run on a Slurm cluster, drop the `execute_cmd` and
`execute_copy_outputs` lines; the defaults submit each step with `srun`.

## Main settings

| Variable | Default | Meaning |
|---|---|---|
| `base_dir`, `project` | — / `ProjectName` | `project_dir` is `{base_dir}/{project}` |
| `fastq_dir`, `meta_dir`, `out_dir` | `{project_dir}/fastq`, `/meta`, `/results` | Input reads, metadata, results |
| `genome_file`, `gff_file` | `{base_dir}/genome.fasta`, `.gff` | Reference genome, annotation |
| `origin_file` | `''` (none) | Origin locations (BED) for velogp, optional |
| `genome_circular` | `True` | Whether chromosomes are circular (`True` or `False`) |
| `sample`, `reference` | — | Sample name and its reference sample; set per target |
| `windowsize` | `10k` | Window size for abundances (`k` = ×1000) |
| `mapping_quality` | `20` | Minimum mapping quality used for abundances and mappability |
| `partition_default` | `compute` | Slurm partition steps are submitted to |
| `cores`, `memory_*_gb`, `time_default` | `32`, …, `1-0` | Resources requested per step |
| `execute_cmd` | `srun -p '{step_partition}' …` | Executor; `'{script_file}'` runs steps locally |
| `execute_copy_outputs` | `True` | Write outputs to `tmp_dir` first, then copy them to their final place |
| `steps` | `genome_index` … `velogp` | Default step list: `genome_index`, `map`, `fastqc`, `mapstats`, `mappability`, `abundance`, `velogp` |

## Steps

Default output names are given relative to `out_dir`. All of them can be
changed through the variable named in each step.

### genome_index

Builds the BWA index of the reference genome. Run it once per genome, as its
own target, before any sample.

* **Input** `genome_file`: reference genome,
  [FASTA](https://en.wikipedia.org/wiki/FASTA_format). Sequence names become the
  chromosome names in all later files.
* **Output** `genome_index_dir` (default `{genome_file}.bwaidx/`): a directory of
  [BWA](https://bio-bwa.sourceforge.net/bwa.shtml) index files with prefix
  `index` (variable `genome_index_prefix`).

### map

Maps read pairs with `bwa mem`, then sorts and indexes the alignments with
`samtools`.

* **Input** `read1_file`, `read2_file` (default `{fastq_dir}/{sample}_r1.fastq.gz`,
  `_r2.fastq.gz`): mates 1 and 2,
  [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format), optionally gzip- or
  bzip2-compressed.
* **Input** `genome_index_dir`: from `genome_index`.
* **Output** `bam_file` (`{sample}.bam`): coordinate-sorted alignments,
  [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf).
* **Output** `bam_index_file` (`{sample}.bai`): BAM index (BAI, defined in the same
  specification).

### fastqc

Runs [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
quality control on the mapped reads.

* **Input** `bam_file`.
* **Output** `fastqc_html_file` (`{sample}_fastqc.html`): FastQC HTML report.
* **Output** `fastqc_zip_file` (`{sample}_fastqc.zip`): ZIP archive with the raw
  FastQC data and plots.

### mapstats

Computes basic mapping statistics, including the median insert size and read
length that `mappability` needs.

* **Input** `bam_file`, `bam_index_file`.
* **Output** `stats_file` (`{sample}.stats`): key/value table, see
  [stats_format.md](stats_format.md).

### mappability

Finds the genomic regions that cannot be mapped reliably with this sample's read
length and insert size. It does this by simulating every possible fragment
and mapping it. Regions can also be excluded by hand.

* **Input** `genome_file`, `genome_index_dir`.
* **Input** `stats_file`: provides read length and insert size.
* **Input** `regions_file` (default `''`, i.e. none): optional
  [BED](https://samtools.github.io/hts-specs/BEDv1.pdf) file with at least 4
  columns. Regions whose name (column 4) is listed in the comma-separated
  `exclude_regions` setting are marked non-mappable.
* **Output** `mappability_file` (`{sample}.mappability.bed`): non-mappable regions,
  see [mappability_format.md](mappability_format.md).

### abundance

Computes coverage in windows of `windowsize` bases for the sample and its
reference, their ratio with standard errors, and per-window mappability
flags.

* **Input** `bam_file`, `bam_index_file`, `mappability_file` of the sample.
* **Input** `reference_bam_file`, `reference_bam_index_file`,
  `reference_mappability_file`: the same files for the sample named in
  `reference`.
* **Output** `abundance_file` (`{sample}.abd.w{windowsize}bp.gz`): see
  [abundance_format.md](abundance_format.md).

To use no reference, add `mfa_no_reference` to `use`, before `mfa_pipeline`
(e.g. `use: [mfa_no_reference, mfa_pipeline]`). It sets the reference files to `''`, so
only the `sample_*` columns are filled. It also switches `velogp` and
`abundance_plot` to those columns.

### abundance_plot

Plots the abundance profile with gnuplot, showing well and badly mappable
windows separately. The plot is set by `abundance_plot_cmd`.

* **Input** `abundance_file`.
* **Output** `abundance_plot` (`{sample}.abd.w{windowsize}bp.pdf`): PDF.

### velocity_oscillation

**Outdated**: use [velogp](#velogp) instead.

Fits a model of replication fork speed to the abundance ratio. The model has
constant or oscillating speed, with or without diffusion. It assumes a
circular genome with a single origin. Genome length and origin come from a
built-in table selected by `genome_name`; currently only `EColiMG1655` is
defined, in `steps/fit_velocity_oscillation.R`.

* **Input** `abundance_file`: uses `midpoint`, `coverage_ratio`,
  `coverage_ratio_stderr`, `ratio_mappable`.
* **Input** the growth rate `growthrate_perhour`, by default looked up in
  `growthrates_file`, see [growthrates_format.md](growthrates_format.md).
* Settings: `genome_name`; `velocity_oscillation_model` (`vt`, `vx` or `vd`;
  default `vt`).
* **Output** `velocity_oscillation_file` (`{sample}.osc.{velocity_oscillation_model}.tab`):
  see [velocity_oscillation_format.md](velocity_oscillation_format.md).

### velogp

Fits a Gaussian-process model of fork speed over time, and optionally of delays at
given genomic features, using [VeloGP](https://github.com/oist/VeloGP). The model and
its options are described in [velogp.md](https://github.com/oist/VeloGP/blob/master/docs/velogp.md).

* **Input** `abundance_file`: the columns used are set by `velogp_abd_columns`.
* **Input** `genome_file` (FASTA).
* **Input** `gff_file`
  ([GFF3](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)
  genome annotation) and `origin_file` (optional, default `''`, i.e. none,
  [BED](https://samtools.github.io/hts-specs/BEDv1.pdf)): the features in which
  origins and vx features are looked up, see
  [features_format.md](https://github.com/oist/VeloGP/blob/master/docs/features_format.md).
* **Input** the growth rate `growthrate_perhour`, by default looked up in
  `growthrates_file`, see [growthrates_format.md](growthrates_format.md).
* Setting `origins` (default `origin_of_replication`): the replication origins, as a
  comma-separated list of locations (`chr:pos`), feature names (from the GFF or BED
  file) or GFF feature types. The default selects all GFF features of type
  `origin_of_replication`.
* Other settings: `velogp_vt_cov`, `velogp_vt_timescale`, `velogp_vt_magnitude`
  (prior for v(t)); `velogp_vx_types`, `velogp_vx_names` (features for which delays
  are estimated); `velogp_errormodel`; `velogp_seed`; `velogp_extra_args`.
  `velogp_suffix` is added to all output names and to the log name.
* **Output** `velogp_abundance_file` (`{sample}.velogp.abd.gz`): observed and fitted
  abundance, see [fitted_abundance_format.md](https://github.com/oist/VeloGP/blob/master/docs/fitted_abundance_format.md).
* **Output** `velogp_vt_file` (`{sample}.velogp.vt.gz`): fitted fork velocity over
  time and replisome trajectories, see
  [trajectories_format.md](https://github.com/oist/VeloGP/blob/master/docs/trajectories_format.md).
* **Output** `velogp_vx_file` (`{sample}.velogp.vx.gz`): feature delays, see
  [vx_feature_delays_format.md](https://github.com/oist/VeloGP/blob/master/docs/vx_feature_delays_format.md).
* **Output** `velogp_abundance_plot`, `velogp_vt_plot` (`.pdf`): plots of the above,
  see [velogp.md](https://github.com/oist/VeloGP/blob/master/docs/velogp.md#plots-and-console-output).
* **Output** `velogp_objects_file` (`{sample}.velogp.pkl`): the fitted model, the
  input for `velomix`, see [objects_format.md](https://github.com/oist/VeloGP/blob/master/docs/objects_format.md).

### velomix

Starts from a fitted VeloGP model and fits the VeloMix mixture model of stochastic
origin firing and replisome stalling at the sites given in `velomix_stall_sites`
(`python -m velogp.velomix`). The model and its options are described in
[velomix.md](https://github.com/oist/VeloGP/blob/master/docs/velomix.md).

* **Input** `velogp_objects_file`: from `velogp`.
* Settings: `velomix_stall_sites`, `velomix_seed`, `velomix_extra_args`.
  `velomix_suffix` is added to all output names.
* **Output** `velomix_probabilities_file` (`{sample}.velomix.p.gz`): mixture
  probabilities, see [probabilities_format.md](https://github.com/oist/VeloGP/blob/master/docs/probabilities_format.md).
* **Output** `velomix_abundance_file`, `velomix_vt_file`, `velomix_vx_file`
  (`.abd.gz`, `.vt.gz`, `.vx.gz`) and plots `velomix_abundance_plot`,
  `velomix_vt_plot` (`.pdf`): same formats as for `velogp`; plots are described in
  [velomix.md](https://github.com/oist/VeloGP/blob/master/docs/velomix.md#plots-and-console-output).
* **Output** `velomix_objects_file` (`{sample}.velomix.pkl`): the fitted model, see
  [objects_format.md](https://github.com/oist/VeloGP/blob/master/docs/objects_format.md).
