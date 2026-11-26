# sra_se_pipeline2

File: `sra_se_pipeline2.yaml`. Paired-end variant:
[sra_pipeline2](sra_pipeline2.md). Alternative that keeps `.sra` files:
[sra_se_pipeline](sra_se_pipeline.md).

## Summary

`sra_se_pipeline2` is the **single-end** variant of [sra_pipeline2](sra_pipeline2.md).
It downloads one or more SRA runs directly with `fasterq-dump` and concatenates
them into a single FASTQ file, `read_file`, as expected by
[mfa_se_pipeline](mfa_se_pipeline.md).

### Use together with mfa_se_pipeline

`sra_se_pipeline2` is not self-contained. It is meant to be listed in `use`
together with [mfa_se_pipeline](mfa_se_pipeline.md), which provides everything it
does not define itself:

* the FASTQ path `read_file` (from `fastq_dir`, `sample`, `base_dir` and
  `project`), which the `fastq` step writes and `map` then reads. Only
  `mfa_se_pipeline` defines `read_file`; `mfa_pipeline` does not;
* the settings every step needs to run: `shell`, `execute_cmd`, `tmp_dir`,
  `log_file`, `step_partition`, `step_cores`, `step_memory_gb` and `step_time`.

`sra_se_pipeline2` and `mfa_se_pipeline` define no settings in common, so their
order in `use` does not matter.

## Quick start

This example downloads a non-replicating reference sample and an exponentially
growing sample made of two runs from SRA. It then analyses them with
`mfa_se_pipeline` up to VeloGP. It assumes this layout:

```
/data/genome.fasta
/data/genome.gff
/data/public/meta/growthrates.tab        # see growthrates_format.md
```

Create `public.yaml`:

```yaml
use:
  - sra_se_pipeline2
  - mfa_se_pipeline
project: public
base_dir: /data                 # project_dir = /data/public
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
      accession: SRR1234567
      steps: [ fastq, map, mapstats, mappability ]
  growing:
    - sample: exp
      accession: SRR1234568,SRR1234569
      reference: stat
      steps: [ fastq, map, mapstats, mappability, abundance, velogp ]
```

velogp looks up the replication origins in `genome.gff`, by default all
features of type `origin_of_replication`. Set `origins` to select them by name
or feature type instead, e.g. `origins: oriC`. If the GFF file does not contain
the origins, list them in a BED file and set `origin_file` to it (see [velogp](mfa_pipeline.md#velogp)).

```bash
arubaito public genome reference growing
```

For each sample, `fastq` creates `/data/public/fastq/<sample>.fastq.gz`; for `exp`
it holds the reads of both runs. The remaining steps are those of
`mfa_se_pipeline`, and their results go to `/data/public/results/`. To only
download, list just `[ fastq ]` in `steps`.

## Steps

### fastq

Runs `fasterq-dump --split-3` for each accession in the job's `tmp_dir`, then
concatenates the single-end reads of all accessions, gzip-compressing them if
the output name ends in `.gz`. Uses one core.

* Setting `accession`: one or more SRA run accessions, comma-separated.
* **Output** `read_file` (`{fastq_dir}/{sample}.fastq.gz`):
  [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format).

`tmp_dir` must have room for the uncompressed reads of all runs in the sample.

Paired reads in the runs are not included in the output; a warning names each
accession that contained them.
