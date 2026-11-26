# Mappability file format

Produced by the `mappability` step of [mfa_pipeline](mfa_pipeline.md) and
[mfa_se_pipeline](mfa_se_pipeline.md) (`steps/mappability.py`). It is read by
the `abundance` step.

The file lists the regions of the genome where reads **cannot** be mapped
reliably. The step simulates every possible fragment of the genome (one per
start position), maps it with `bwa mem` using the sample's read length and
fragment length, and marks a region as non-mappable if it contains:

* a fragment that was mapped to the wrong place, with mapping quality below
  `mapping_quality`, or not as a full-length match; or
* the place such a fragment was wrongly mapped to; or
* a region from `regions_file` whose name is listed in `exclude_regions`.

Because the result depends on read and fragment length, each sample gets its
own mappability file.

## Format

A [BED](https://samtools.github.io/hts-specs/BEDv1.pdf) file with three columns
(`chrom`, `start`, `end`), tab-separated, using 0-based half-open coordinates.
Each line is one **non-mappable** region. Regions are sorted within each
chromosome.

The file starts with `#` comment lines that record how it was made. Readers
skip them:

```
# Command used to create this file:
#   .../steps/mappability.py --genome genome.fasta --genome-index genome.fasta.bwaidx/index --circular --paired --fragment-length 312 --read-length 150 ...
# Chromosomes:
#   NC_000913.3: 4641652
NC_000913.3	15387	16801
NC_000913.3	223408	228875
```

The `# Chromosomes:` block lists every sequence in the genome with its length.
