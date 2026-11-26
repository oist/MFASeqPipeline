# sra_pipeline2

File: `sra_pipeline2.yaml`. Single-end variant:
[sra_se_pipeline2](sra_se_pipeline2.md). Alternative that keeps `.sra` files:
[sra_pipeline](sra_pipeline.md).

## Summary

`sra_pipeline2` downloads **paired-end** runs from the NCBI
[Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra/docs/) straight to
FASTQ with `fasterq-dump` (from [sra-tools](https://github.com/ncbi/sra-tools)),
in a single `fastq` step. It does not keep `.sra` files. A sample may consist of
several runs: `accession` may be a comma-separated list. The runs are
downloaded in parallel and concatenated into one pair of FASTQ files.

### Use together with mfa_pipeline

`sra_pipeline2` is not self-contained. It is meant to be listed in `use` together
with [mfa_pipeline](mfa_pipeline.md), which provides everything it does not
define itself:

* the FASTQ paths `read1_file` and `read2_file` (from `fastq_dir`, `sample`,
  `base_dir` and `project`), which the `fastq` step writes and `map` then reads;
* the settings every step needs to run: `shell`, `execute_cmd`, `tmp_dir`,
  `log_file`, `step_partition`, `step_cores`, `step_memory_gb` and `step_time`.

`sra_pipeline2` and `mfa_pipeline` define no settings in common, so their order
in `use` does not matter.

## Quick start

This example downloads a non-replicating reference sample and an exponentially
growing sample made of two runs from SRA. It then analyses them with
`mfa_pipeline` up to VeloGP. It assumes this layout:

```
/data/genome.fasta
/data/genome.gff
/data/public/meta/growthrates.tab        # see growthrates_format.md
```

Create `public.yaml`:

```yaml
use:
  - sra_pipeline2
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

For each sample, `fastq` creates `/data/public/fastq/<sample>_r1.fastq.gz` and
`_r2.fastq.gz`; for `exp` they hold the reads of both runs. The remaining steps
are those of `mfa_pipeline`, and their results go to `/data/public/results/`. To
only download, list just `[ fastq ]` in `steps`.

## Steps

### fastq

Runs `fasterq-dump --split-3` for each accession in the job's `tmp_dir`. It then
concatenates mates 1 and mates 2 across accessions, gzip-compressing them if the
output name ends in `.gz`. Unpaired reads are dropped with a warning. Requests
two cores.

* Setting `accession`: one or more SRA run accessions, comma-separated.
* **Output** `read1_file`, `read2_file` (`{fastq_dir}/{sample}_r1.fastq.gz`,
  `_r2.fastq.gz`): mates 1 and 2,
  [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format).

`tmp_dir` must have room for the uncompressed reads of all runs in the sample.
