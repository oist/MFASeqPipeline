# sra_pipeline

File: `sra_pipeline.yaml`. Single-end variant:
[sra_se_pipeline](sra_se_pipeline.md). Alternative without `.sra` files, which
also supports several runs per sample: [sra_pipeline2](sra_pipeline2.md).

## Summary

`sra_pipeline` downloads **paired-end** sequencing runs from the NCBI
[Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra/docs/) and converts
them to the gzipped FASTQ files that [mfa_pipeline](mfa_pipeline.md) expects
(`read1_file`, `read2_file`). It works in two steps: `sra` fetches the
`.sra` archive with `prefetch`, and `fastq` extracts it with `fastq-dump`. Both tools
come from [sra-tools](https://github.com/ncbi/sra-tools). Each sample corresponds
to one run accession (variable `accession`).

### Use together with mfa_pipeline

`sra_pipeline` is not self-contained. It is meant to be listed in `use` together
with [mfa_pipeline](mfa_pipeline.md), which provides everything it does not
define itself:

* `project_dir` (from `base_dir` and `project`), which sets the location of the
  `.sra` files (`sra_dir` = `{project_dir}/sra`);
* the FASTQ paths `read1_file` and `read2_file` (from `fastq_dir` and `sample`),
  which the `fastq` step writes and `map` then reads;
* the settings every step needs to run: `shell`, `execute_cmd`, `tmp_dir`,
  `log_file`, `step_partition`, `step_cores`, `step_memory_gb` and `step_time`.

`sra_pipeline` and `mfa_pipeline` define no settings in common, so their order
in `use` does not matter.

## Quick start

This example downloads a non-replicating reference run and an exponentially
growing run from SRA, and analyses them with `mfa_pipeline` up to VeloGP.
It assumes this layout:

```
/data/genome.fasta
/data/genome.gff
/data/public/meta/growthrates.tab        # see growthrates_format.md
```

Create `public.yaml`:

```yaml
use:
  - sra_pipeline
  - mfa_pipeline
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
      steps: [ sra, fastq, map, mapstats, mappability ]
  growing:
    - sample: exp
      accession: SRR1234568
      reference: stat
      steps: [ sra, fastq, map, mapstats, mappability, abundance, velogp ]
```

velogp looks up the replication origins in `genome.gff`, by default all
features of type `origin_of_replication`. Set `origins` to select them by name
or feature type instead, e.g. `origins: oriC`. If the GFF file does not contain
the origins, list them in a BED file and set `origin_file` to it (see [velogp](mfa_pipeline.md#velogp)).

```bash
arubaito public genome reference growing
```

For each sample, `sra` and `fastq` create `/data/public/sra/<accession>.sra` and
`/data/public/fastq/<sample>_r1.fastq.gz` / `_r2.fastq.gz`. The remaining steps
are those of `mfa_pipeline`, and their results go to `/data/public/results/`. To
only download, list just `[ sra, fastq ]` in `steps`.

## Steps

Both steps use one core.

### sra

Downloads the run with `prefetch` and verifies it.

* Setting `accession`: SRA run accession (`SRR…`, `ERR…` or `DRR…`).
* **Output** `sra_file` (`{sra_dir}/{accession}.sra`, with `sra_dir` =
  `{project_dir}/sra`): SRA archive in NCBI's binary SRA format, readable with
  sra-tools.

### fastq

Extracts the reads with `fastq-dump --split-3` and compresses them with gzip if
the output name ends in `.gz`. Unpaired reads that `--split-3` puts into a
separate file are dropped, with a warning.

* **Input** `sra_file`: from `sra`.
* **Output** `read1_file`, `read2_file` (`{fastq_dir}/{sample}_r1.fastq.gz`,
  `_r2.fastq.gz`): mates 1 and 2,
  [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format).
