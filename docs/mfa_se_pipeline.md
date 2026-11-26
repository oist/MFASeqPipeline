# mfa_se_pipeline

File: `mfa_se_pipeline.yaml` (`use: mfa_pipeline`). Paired-end variant:
[mfa_pipeline](mfa_pipeline.md).

## Summary

`mfa_se_pipeline` is the **single-end** variant of [mfa_pipeline](mfa_pipeline.md).
It sets `paired_end: False`, reads one FASTQ file per sample, and replaces the `map`
and `mappability` steps. Everything else (settings, other steps, output formats)
comes from `mfa_pipeline`, so that documentation applies unchanged. With
single-end reads, the `abundance` step counts only the bases covered by each
read, since the full fragment is unknown.

## Quick start

Assume this layout (the same as in the
[mfa_pipeline quick start](mfa_pipeline.md#quick-start), but with one FASTQ file
per sample):

```
/data/genome.fasta
/data/genome.gff
/data/myproject/fastq/stat.fastq.gz
/data/myproject/fastq/exp.fastq.gz
/data/myproject/meta/growthrates.tab     # see growthrates_format.md
```

Create `myproject.yaml`:

```yaml
use: mfa_se_pipeline
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
      steps: [ map, mapstats, mappability, abundance, velogp ]
```

velogp looks up the replication origins in `genome.gff`, by default all
features of type `origin_of_replication`. Set `origins` to select them by name
or feature type instead, e.g. `origins: oriC`. If the GFF file does not contain
the origins, list them in a BED file and set `origin_file` to it (see [velogp](mfa_pipeline.md#velogp)).

```bash
arubaito myproject genome reference growing
```

Results go to `/data/myproject/results/`. For example, `exp.abd.w10kbp.gz` holds
the abundances and `exp.velogp.vt.gz` the fitted fork speed over time.

## Steps

Only the steps that differ from [mfa_pipeline](mfa_pipeline.md#steps) are listed.

### map

Maps single-end reads with `bwa mem`, then sorts and indexes the alignments.

* **Input** `read_file` (default `{fastq_dir}/{sample}.fastq.gz`):
  [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format), optionally gzip- or
  bzip2-compressed.
* **Input** `genome_index_dir`: from `genome_index`. (Unlike in `mfa_pipeline`,
  it is not listed as a formal input. It is only read through
  `genome_index_prefix`, so it is not staged with `execute_copy_inputs`.)
* **Output** `bam_file` (`{sample}.bam`): sorted
  [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf).
* **Output** `bam_index_file` (`{sample}.bai`): BAM index.

### mappability

Same as in `mfa_pipeline`, but the simulated fragments are single reads of the
median read length from `stats_file`. Each read is mapped once, and the
fragment length plays no role.

* **Input** `genome_file`, `genome_index_dir`, `regions_file`, `stats_file`.
* **Output** `mappability_file` (`{sample}.mappability.bed`): see
  [mappability_format.md](mappability_format.md).
