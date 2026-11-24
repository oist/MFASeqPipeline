# sra_se_pipeline

File: `sra_se_pipeline.yaml` (`use: sra_pipeline`). Paired-end variant:
[sra_pipeline](sra_pipeline.md). Alternative without `.sra` files, which also
supports several runs per sample: [sra_se_pipeline2](sra_se_pipeline2.md).

## Summary

`sra_se_pipeline` is the **single-end** variant of [sra_pipeline](sra_pipeline.md).
The `sra` step is the same. The `fastq` step is replaced so that it produces one
FASTQ file, `read_file`, which is what [mfa_se_pipeline](mfa_se_pipeline.md)
expects.

### Use together with mfa_se_pipeline

`sra_se_pipeline` is not self-contained. It is meant to be listed in `use`
together with [mfa_se_pipeline](mfa_se_pipeline.md), which provides everything it
does not define itself:

* `project_dir` (from `base_dir` and `project`), which sets the location of the
  `.sra` files (`sra_dir` = `{project_dir}/sra`);
* the FASTQ path `read_file` (from `fastq_dir` and `sample`), which the `fastq`
  step writes and `map` then reads. Only `mfa_se_pipeline` defines `read_file`;
  `mfa_pipeline` does not;
* the settings every step needs to run: `shell`, `execute_cmd`, `tmp_dir`,
  `log_file`, `step_cores`, `step_memory_gb` and `step_time`.

`sra_se_pipeline` and `mfa_se_pipeline` define no settings in common, so their
order in `use` does not matter.

## Quick start

This example downloads a non-replicating reference run and an exponentially
growing run from SRA, and analyses them with `mfa_se_pipeline` up to VeloGP.
It assumes this layout:

```
/data/genome.fasta
/data/genome.gff
/data/origins.bed                        # origins, for velogp
/data/public/meta/growthrates.tab        # see growthrates_format.md
```

Create `public.yaml`:

```yaml
use: [ sra_se_pipeline, mfa_se_pipeline ]
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

```bash
./arubaito.py public genome reference growing
```

For each sample, `sra` and `fastq` create `/data/public/sra/<accession>.sra` and
`/data/public/fastq/<sample>.fastq.gz`. The remaining steps are those of
`mfa_se_pipeline`, and their results go to `/data/public/results/`. To only
download, list just `[ sra, fastq ]` in `steps`.

## Steps

Both steps use one core.

### sra

As in [sra_pipeline](sra_pipeline.md#sra): downloads `{sra_dir}/{accession}.sra`.

### fastq

Extracts the reads with `fastq-dump --split-3` and compresses them with gzip if
the output name ends in `.gz`. If the run contains paired reads, they are dropped
with a warning.

* **Input** `sra_file`: from `sra`.
* **Output** `read_file` (`{fastq_dir}/{sample}.fastq.gz`):
  [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format).
