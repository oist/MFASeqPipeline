# rnaseq_pipeline

File: `rnaseq_pipeline.yaml`.

## Summary

`rnaseq_pipeline` adds a single step, `genecounts`, that counts RNA-seq fragments
per gene with [featureCounts](https://subread.sourceforge.net/featureCounts.html)
(Subread package). Counting assumes **paired-end, reverse-stranded** libraries
(e.g. dUTP/TruSeq stranded).

### Use together with mfa_pipeline

`rnaseq_pipeline` is not self-contained: it has no mapping steps and defines
only `genecounts_file` and the `genecounts` step. It is meant to be listed in
`use` together with [mfa_pipeline](mfa_pipeline.md), which provides:

* the steps that produce its inputs: `genome_index`, `map` (reads → `bam_file`)
  and `mapstats` (→ `stats_file`);
* the file paths `bam_file`, `stats_file`, `genome_file`, `gff_file` and
  `out_dir` (from `base_dir` and `project`);
* the settings every step needs to run: `shell`, `execute_cmd`, `tmp_dir`,
  `log_file`, `step_cores`, `step_memory_gb` and `step_time`.

`rnaseq_pipeline` and `mfa_pipeline` define no settings in common, so their order
in `use` does not matter. The MFA-specific steps of `mfa_pipeline` (`mappability`,
`abundance`, …) are simply not listed in `steps`.

## Quick start

This example maps one RNA-seq sample and counts fragments per gene. It assumes
the genome and its annotation are at `/data/genome.fasta` and `/data/genome.gff`
(the `mfa_pipeline` defaults), and the reads are in
`/data/myproject/fastq/rna1_r1.fastq.gz` / `rna1_r2.fastq.gz`. Create
`rnaseq.yaml`:

```yaml
use:
  - rnaseq_pipeline
  - mfa_pipeline
project: myproject
base_dir: /data                 # project_dir = /data/myproject
# Run steps locally instead of with Slurm
execute_cmd: '{script_file}'
execute_copy_outputs: False
cores: 4
targets:
  genome:
    - target_label: genome
      steps: [ genome_index ]
  samples:
    - sample: rna1
      steps: [ map, mapstats, genecounts ]
```

```bash
arubaito rnaseq genome samples
```

The counts are written to `/data/myproject/results/rna1.genecounts.tab`.

To download the genome from NCBI instead, add
[refseq_pipeline](refseq_pipeline.md) **before** `mfa_pipeline` in `use`. The
earlier entry wins, so its `genome_file` and `gff_file` replace the defaults:

```yaml
use:
  - refseq_pipeline
  - rnaseq_pipeline
  - mfa_pipeline
refseq: GCF_000005845.2
targets:
  genome:
    - target_label: genome
      steps: [ genome_download, genome_index ]
  # ... samples as above
```

## Steps

### genecounts

Counts fragments overlapping each `gene` feature of the annotation.

* **Input** `bam_file` (`{out_dir}/{sample}.bam`): sorted paired-end alignments,
  [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf), e.g. from `map`.
* **Input** `stats_file` (`{out_dir}/{sample}.stats`): from `mapstats`, see
  [stats_format.md](stats_format.md). The median `insert_size` sets the accepted
  fragment length range: insert size / 5 to insert size × 5.
* **Input** `gff_file`: annotation,
  [GFF3](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md).
  Features of type `gene` are counted and identified by their `ID` attribute.
* **Input** `genome_file`: listed for completeness but not used by the script.
* **Output** `genecounts_file` (`{out_dir}/{sample}.genecounts.tab`): featureCounts
  count table. It is tab-separated, starts with a `#` comment line holding the
  program call, and has a header line with the columns `Geneid`, `Chr`, `Start`,
  `End`, `Strand`, `Length`, then the extra attributes `gene`, `gene_synonym`,
  `Dbxref`, `gene_biotype`, and a last column with the count (named after the
  BAM file). featureCounts also writes `{genecounts_file}.summary`, with the
  numbers of assigned and unassigned fragments by reason. The pipeline does not
  track this file.

Counting rules (featureCounts options): pairs are counted as fragments
(`-p --countReadPairs`), reverse-stranded (`-s 2`), with mapping quality ≥ 10
(`-Q 10`). Both mates must be mapped (`-B`), chimeric fragments are not counted
(`-C`), and the fragment length must be within the range above (`-P -d -D`). A
fragment overlapping several genes is counted for each of them (`-O`).
