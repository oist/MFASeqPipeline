# Mapping statistics file format

Produced by the `mapstats` step of [mfa_pipeline](mfa_pipeline.md) and
[mfa_se_pipeline](mfa_se_pipeline.md) (`steps/mapping_statistics.py`). The
`mappability` and `genecounts` steps read the read length and insert size from
it.

Tab-separated text with one `key<TAB>value` pair per line and no header.
Only primary alignments are counted, so each read is counted once.

| Key | Description |
|---|---|
| `total_reads` | Number of reads (mates count separately) |
| `mapped_reads` | Number of mapped reads |
| `properly_paired_reads` | Number of reads flagged as part of a proper pair (0 for single-end data) |
| `insert_size` | Median absolute template length (fragment length) of properly paired reads; `nan` if there are none |
| `read_length` | Median read length |

Example:

```
total_reads	9123456
mapped_reads	9034567
properly_paired_reads	8912345
insert_size	312
read_length	150
```

Other steps extract values like this: `grep insert_size "{stats_file}" | cut -f 2`.
